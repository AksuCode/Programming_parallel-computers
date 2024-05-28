#include "llm.h"
#include "ppc.h"
#include <fstream>
#include <utility>

namespace ppc {
std::vector<float> llama_reference(LLamaConfig config, LLamaParameters params, std::vector<token_t> prompt);

std::pair<LLamaConfig, LLamaParameters> make_test_model(LLamaConfig config,
                                                        int multiplier,
                                                        const std::string &modifier,
                                                        ppc::random rng,
                                                        bool for_test);

bool test_against_ref_implementation(std::ostream *stream, LLamaConfig base_config, ppc::random rng,
                                     int multiplier, const std::string &model_modifier,
                                     std::vector<token_t> input, const std::vector<float> &logits);
} // namespace ppc

int main(int argc, const char **argv) {
    const char *ppc_output = std::getenv("PPC_OUTPUT");
    int ppc_output_fd = 0;
    if (ppc_output) {
        ppc_output_fd = std::stoi(ppc_output);
    }
    if (ppc_output_fd <= 0) {
        ppc_output_fd = 1;
    }

    std::unique_ptr<ppc::fdostream> stream = std::unique_ptr<ppc::fdostream>(new ppc::fdostream(ppc_output_fd));

    argc--;
    argv++;
    if (argc < 1 || argc > 2) {
        std::cerr << "Invalid usage" << std::endl;
        return 1;
    }

    bool test = false;
    if (argv[0] == std::string("--test")) {
        test = true;
        argc--;
        argv++;
    }

    std::ifstream input_file(argv[0]);
    if (!input_file) {
        std::cerr << "Failed to open input file" << std::endl;
        return 2;
    }

    std::string time_out;
    CHECK_READ(input_file >> time_out);
    if (time_out == "timeout") {
        input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::string model_modifier;
    int model_multiplier;
    int num_tokens, seed;
    auto read_cfg = [&](const std::string &expect, auto &target) {
        std::string line;
        std::getline(input_file, line);
        std::stringstream line_data(line);
        CHECK_READ(line_data >> line);
        if (line != expect) {
            std::cerr << "Invalid test spec: Expected " << expect << ". got " << line << std::endl;
            std::exit(1);
        }
        // no CHECK_READ macro; the check below does the same thing, but gives a better error message
        line_data >> target;
        if (!line_data) {
            std::cerr << "Error while reading test spec " << expect << "." << std::endl;
            std::exit(1);
        }
    };

    LLamaConfig base_config;
    read_cfg("seed", seed);
    read_cfg("dim", base_config.dim);
    read_cfg("n_layers", base_config.n_layers);
    read_cfg("hidden_dim", base_config.hidden_dim);
    read_cfg("n_heads", base_config.n_heads);
    read_cfg("vocab_size", base_config.vocab_size);
    read_cfg("seq_len", base_config.seq_len);
    read_cfg("multiplier", model_multiplier);
    read_cfg("modifier", model_modifier);
    read_cfg("num_tokens", num_tokens);

    // check that we can actually use the requested multiplier
    if (base_config.n_layers % model_multiplier != 0) {
        std::cerr << "Error, invalid multiplier " << model_multiplier << "  for model with " << base_config.n_layers << "layers" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // validate that we keep our promises
    if (base_config.hidden_dim % 16 != 0) {
        std::cerr << "Error, invalid hidden_dim " << base_config.hidden_dim << ": not a multiple of 16." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (base_config.dim % 16 != 0) {
        std::cerr << "Error, invalid dim " << base_config.dim << ": not a multiple of 16." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (base_config.head_size() % 16 != 0) {
        std::cerr << "Error, invalid head size " << base_config.head_size() << " ("
                  << base_config.dim << " / " << base_config.n_heads << "): not a multiple of 16." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::vector<token_t> prompt_tokens(num_tokens);
    ppc::random rng(seed);
    for (int i = 0; i < num_tokens; ++i) {
        prompt_tokens[i] = (token_t)rng.get_int64(0, base_config.vocab_size);
    }

    auto [config, params] = ppc::make_test_model(base_config, model_multiplier, model_modifier, rng, false);

    std::vector<float> logits(prompt_tokens.size() * config.vocab_size);

    ppc::setup_cuda_device();
    ppc::perf timer;
    timer.start();

    // first, init the model
    llm(config, std::move(params), prompt_tokens, logits);

    // then do the actual run
    timer.stop();
    timer.print_to(*stream);
    ppc::reset_cuda_device();

    if (test) {
        bool pass = ppc::test_against_ref_implementation(stream.get(), base_config, rng,
                                                         model_multiplier, model_modifier,
                                                         prompt_tokens, logits);
        if (!pass) {
            *stream << "seed\t";
            *stream << seed << "\n";
            *stream << "length\t";
            *stream << num_tokens << "\n";

            // give some context rgd. the model
            *stream << "config.vocab_size\t" << config.vocab_size << "\n";
            *stream << "config.dim\t" << config.dim << "\n";
            *stream << "config.hidden_dim\t" << config.hidden_dim << "\n";
            *stream << "config.n_layers\t" << config.n_layers << "\n";
            *stream << "config.seq_len\t" << config.seq_len << "\n";
            *stream << "config.n_heads\t" << config.n_heads << "\n";

            *stream << '\n';
            *stream << std::endl;
            exit(0);
        }

        *stream << "result\tpass\n";
    } else {
        *stream << "result\tdone\n";
    }
    *stream << std::endl;
}

std::vector<float> rand_vec(ppc::random &rng, std::size_t size) {
    std::vector<float> vec(size);
    std::generate(vec.begin(), vec.end(), [&]() {
        return rng.get_float(-1, 1);
    });
    return vec;
}

std::pair<LLamaConfig, LLamaParameters> ppc::make_test_model(
    LLamaConfig config,
    int multiplier,
    const std::string &modifier,
    ppc::random rng,
    bool for_test) {
    LLamaParameters weights;
    weights.TokenEmbeddingMatrix = rand_vec(rng, config.vocab_size * config.dim);

    // under-provision the LayerWeights; only config.n_layers / multiplier will be non-trivial, to speed up testing
    weights.LayerWeights.resize(config.n_layers / multiplier);

    auto make_data = [&](auto &&ptr_to_member, int num) {
        for (auto &layer : weights.LayerWeights) {
            layer.*ptr_to_member = rand_vec(rng, num);
        }
    };

    int head_size = config.dim / config.n_heads;
    make_data(&LLamaLayer::rms_attention, config.dim);
    make_data(&LLamaLayer::query_weight_matrix, config.dim * (config.n_heads * head_size));
    make_data(&LLamaLayer::key_weight_matrix, config.dim * (config.n_heads * head_size));
    make_data(&LLamaLayer::value_weight_matrix, config.dim * (config.n_heads * head_size));
    make_data(&LLamaLayer::out_weight_matrix, (config.n_heads * head_size) * config.dim);

    make_data(&LLamaLayer::rms_feed_forward, config.dim);
    make_data(&LLamaLayer::feed_forward_w1, config.dim * config.hidden_dim);
    make_data(&LLamaLayer::feed_forward_w2, config.hidden_dim * config.dim);
    make_data(&LLamaLayer::feed_forward_w3, config.dim * config.hidden_dim);

    weights.RmsFinal = rand_vec(rng, config.dim);
    // unused
    rand_vec(rng, config.seq_len * head_size);
    weights.TokenOutputMatrix = rand_vec(rng, config.vocab_size * config.dim);

    if (modifier == "default") {
    } else if (modifier == "disable-attention") {
        // set out_weight_matrix to zero, effectively disabling the attention module
        for (auto &l : weights.LayerWeights) {
            l.out_weight_matrix.assign(l.out_weight_matrix.size(), 0.f);
        }
    } else if (modifier == "disable-feed-forward") {
        // set ffn w2 to zero, effectively disabling the feed-forward module
        for (auto &l : weights.LayerWeights) {
            l.feed_forward_w2.assign(l.feed_forward_w2.size(), 0.f);
        }
    } else {
        std::cerr << "Invalid modifier '" << modifier << "'" << std::endl;
        std::exit(1);
    }

    // if we have a reference model, the effective number of layers is divided by the multiplier
    // if we have the student's model, add in some dummy layers that effectively get remove by
    // having all their outputs multiplied by zero.
    if (!for_test) {
        for (int k = 0; k < multiplier; ++k) {
            for (int i = 0; i < config.n_layers / multiplier; ++i) {
                weights.LayerWeights.push_back(weights.LayerWeights[i]);
                auto &l = weights.LayerWeights.back();
                l.out_weight_matrix.assign(l.out_weight_matrix.size(), 0.f);
                l.feed_forward_w2.assign(l.feed_forward_w2.size(), 0.f);
            }
        }
    } else {
        config.n_layers = config.n_layers / multiplier;
    }

    return std::make_pair(config, std::move(weights));
}

bool ppc::test_against_ref_implementation(std::ostream *stream, LLamaConfig base_config, ppc::random rng, int model_multiplier,
                                          const std::string &model_modifier,
                                          std::vector<token_t> input, const std::vector<float> &logits) {
    // Load and run the reference model
    auto [config, params] = ppc::make_test_model(base_config, model_multiplier, model_modifier, rng, true);
    auto ref_result = ppc::llama_reference(config, std::move(params), std::move(input));

    const float threshold = 1e-2;
    float max_error = threshold;
    unsigned max_error_pos = -1;
    for (unsigned i = 0; i < logits.size(); ++i) {
        if (std::abs(ref_result[i] - logits[i]) > max_error) {
            max_error = std::abs(ref_result[i] - logits[i]);
            max_error_pos = i;
        }
    }

    if (max_error > threshold) {
        *stream << "result\tfail\n";
        *stream << "location_pos\t" << max_error_pos / config.vocab_size << "\n";
        *stream << "location_tok\t" << max_error_pos % config.vocab_size << "\n";
        *stream << "threshold\t" << threshold << "\n";
        *stream << "max_error\t" << max_error << "\n";
        return false;
    }

    return true;
}
