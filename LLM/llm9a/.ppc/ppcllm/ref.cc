#include "llm.h"
#include "omp.h"
#include <numeric>
#include <vector>

using namespace utils;

namespace {

inline void matmul(float *out, const float *x, const float *w, int n, int d) {
// when running grading on student's computers, we're working on a system that is in interactive use.
// trying to get _all_ threads here results in terrible performance. Since this loop is very short,
// if just one of the threads gets stolen away by the OS for another task, everything has to wait.
// Therefore, we leave one thread unused.
#pragma omp parallel for num_threads(std::max(1, omp_get_num_threads() - 1))
    for (int i = 0; i < d; i++) {
        // transform-reduce is a slightly faster version of std::inner_product
        // it is not allowed in student code
        out[i] = std::transform_reduce(x, x + n, w + i * n, 0.0);
    }
}

class LLamaModel {
  public:
    LLamaModel(LLamaConfig config, LLamaParameters params);
    std::vector<float> predict(int position, token_t token);

  private:
    // helpers running part of the model
    void multi_head_attention(const LLamaLayer &layer,
                              float *activation, float *attention,
                              int position, int layer_id,
                              float *query);

    LLamaConfig Config;
    LLamaParameters Parameters;

    std::vector<float> KeyCache;   // (layer, seq_len, dim)
    std::vector<float> ValueCache; // (layer, seq_len, dim)
};

LLamaModel::LLamaModel(LLamaConfig config, LLamaParameters params) : Config(config), Parameters(std::move(params)) {
    KeyCache.resize(Config.n_layers * Config.seq_len * Config.dim);
    ValueCache.resize(Config.n_layers * Config.seq_len * Config.dim);
}

std::vector<float> LLamaModel::predict(int position, token_t token) {
    // a few convenience variables
    int dim = Config.dim;
    int hidden_dim = Config.hidden_dim;

    // initialize the activation by the embedding of the current token
    std::vector<float> activation(Parameters.TokenEmbeddingMatrix.begin() + (int)token * dim,
                                  Parameters.TokenEmbeddingMatrix.begin() + ((int)token + 1) * dim);

    // scratch buffers to use during the computations
    std::vector<float> buffer(dim);
    std::vector<float> buffer2(dim);
    std::vector<float> hidden_buffer(hidden_dim);
    std::vector<float> hidden_buffer2(hidden_dim);
    std::vector<float> query(dim);
    std::vector<float> attention(Config.n_heads * Config.seq_len);

    // forward all the layers
    for (int l = 0; l < Config.n_layers; l++) {
        auto &layer = Parameters.LayerWeights[l];

        // attention rmsnorm
        rmsnorm(buffer.data(),
                activation.data(),
                layer.rms_attention.data(),
                dim);

        multi_head_attention(layer, buffer.data(), attention.data(), position, l, query.data());

        // final matmul to get the output of the attention
        matmul(buffer2.data(), buffer.data(), layer.out_weight_matrix.data(), dim, dim);

        // residual connection back into x
        for (int i = 0; i < dim; i++) {
            activation[i] += buffer2[i];
        }

        rmsnorm(buffer.data(), activation.data(), layer.rms_feed_forward.data(), dim);

        matmul(hidden_buffer.data(), buffer.data(), layer.feed_forward_w1.data(), dim, hidden_dim);
        matmul(hidden_buffer2.data(), buffer.data(), layer.feed_forward_w3.data(), dim, hidden_dim);
        swiglu(hidden_buffer.data(), hidden_buffer.data(), hidden_buffer2.data(), hidden_dim);
        // final matmul to get the output of the ffn
        matmul(buffer.data(), hidden_buffer.data(), layer.feed_forward_w2.data(), hidden_dim, dim);

        // residual connection back into x
        for (int i = 0; i < dim; i++) {
            activation[i] += buffer[i];
        }
    }

    // final rmsnorm
    rmsnorm(activation.data(), activation.data(), Parameters.RmsFinal.data(), dim);

    // classifier into logits
    std::vector<float> logits(Config.vocab_size);
    matmul(logits.data(), activation.data(), Parameters.TokenOutputMatrix.data(), dim, Config.vocab_size);
    return logits;
}

void LLamaModel::multi_head_attention(
    const LLamaLayer &layer,
    float *activation, float *attention,
    int position, int layer_id,
    float *query) {
    int dim = Config.dim;
    int head_size = Config.head_size();

    // key and value point to the kv cache
    int offset = layer_id * Config.seq_len * dim; // kv cache layer offset for convenience

    float *key = KeyCache.data() + offset + position * dim;
    float *value = ValueCache.data() + offset + position * dim;

    // qkv matmuls for this position. key and value are generated directly at the desired position
    // inside the cache
    matmul(query,
           activation,
           layer.query_weight_matrix.data(),
           dim, dim);
    matmul(key,
           activation,
           layer.key_weight_matrix.data(),
           dim, dim);
    matmul(value,
           activation,
           layer.value_weight_matrix.data(),
           dim, dim);

    rope(Config, query, key, position);

    const float *key_base = KeyCache.data() + offset;
    const float *value_base = ValueCache.data() + offset;

    // multihead attention. iterate over all heads
#pragma omp parallel for
    for (int h = 0; h < Config.n_heads; h++) {
        // get the query vector for this head
        float *q = query + h * head_size;
        // attention scores for this head
        float *att = attention + h * Config.seq_len;

        calculate_attention(Config, att, q, position,
                            key_base + h * head_size);
        lookup_with_attention(Config, att,
                              activation + h * head_size,
                              position,
                              value_base + h * head_size);
    }
}

} // namespace

namespace ppc {
std::vector<float> llama_reference(LLamaConfig config, LLamaParameters params, std::vector<token_t> prompt);
} // namespace ppc

std::vector<float> ppc::llama_reference(LLamaConfig config, LLamaParameters params, std::vector<token_t> prompt) {
    LLamaModel model(config, std::move(params));
    std::vector<float> result;
    result.reserve(config.seq_len * prompt.size());
    for (unsigned i = 0; i < prompt.size(); ++i) {
        auto pred = model.predict(i, prompt[i]);
        result.insert(result.end(), pred.begin(), pred.end());
    }
    return result;
}
