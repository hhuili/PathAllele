// Emit unpadded edits from vcflib's original BiWFA engine, before VCF serialization.
// One complete REF/ALT pair per input row. GTs are deliberately not interpreted.
#include "Variant.h"
#include "vcf-wfa.h"
#include <iostream>
using namespace vcflib;
int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "Usage: align_alleles input.vcf\n"; return 2; }
    VariantCallFile input;
    input.open(argv[1]);
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
}
