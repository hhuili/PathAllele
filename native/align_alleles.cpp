// Emit unpadded edits from vcflib's original WFA2 backend, before VCF serialization.
// One complete REF/ALT pair per input row. GTs are deliberately not interpreted.
#include "Variant.h"
#include "vcf-wfa.h"
#include <iostream>
using namespace vcflib;
static int align_file(const char* path) {
    VariantCallFile input;
    input.open(path);
    if (!input.is_open()) return 2;
    auto params = wavefront_aligner_attr_default;
    params.distance_metric = gap_affine_2p;
    params.alignment_scope = compute_alignment;
    params.affine2p_penalties = {0, 19, 39, 3, 81, 1};
    WfaVariant row(input);
    while (input.getNextVariant(row)) {
        if (row.alt.size() != 1) { std::cerr << "Expected one ALT\n"; return 2; }
        auto result = row.wfa_parsedAlternates(false, false, false, "", "", &params, 17, 1000, 1, false);
        const auto& aligned = result.at(row.alt.front());
        for (const auto& edit : aligned.first) {
            if (edit.ref == edit.alt) continue;
            std::cout << row.id << '\t' << edit.position << '\t'
                      << (edit.ref.empty() ? "-" : edit.ref) << '\t'
                      << (edit.alt.empty() ? "-" : edit.alt) << '\t'
                      << aligned.second << '\n';
        }
    }
    return 0;
}

// Never allow C++ exceptions to unwind across the C ABI.
extern "C" int pathallele_align_file(const char* path) {
    try { return align_file(path); }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
    catch (...) { std::cerr << "Unknown native alignment error\n"; return 2; }
}
