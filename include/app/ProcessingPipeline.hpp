#ifndef JGAP_PROCESSINGPIPELINE_HPP
#define JGAP_PROCESSINGPIPELINE_HPP
#include <utility>

#include "data/DataNode.hpp"
#include "instructions/PipelineInstruction.hpp"

namespace jgap {
    class ProcessingPipeline {
    public:
        static ProcessingPipeline fromDataNode(const DataNode &node);

        ProcessingPipeline(std::vector<std::shared_ptr<PipelineInstruction>> instructions)
            : instructions_(std::move(instructions)) {}
        void execute();
    private:
        std::vector<std::shared_ptr<PipelineInstruction>> instructions_;
    };
}

/*
- load_xyz:
    file: abc.xyz (=> label = abc)
(
- load_xyz:
    file: abc.xyz => ERROR label taken
)
- modify_data:
    label: abc.xyz
    mod: [delete_virials, delete_energies]
- fit:
    type: isolated_atom
    label: ia1 (-> undef => label=isolated_atom)
    run: true
    => current pot
- potential:
    (label: zbl)
    type: zbl
    => current pot = [ia1, zbl] ?
- potential:
    from_file: pot1.yaml
    => current pot = [ia1, zbl, pot1.yaml]  is it a good idea ?
- fit:
    (fit_?)label: gap2b
    type: gap
    regularization:
        type: multiplier_per_config_type
        E: 0.1
        is:
            isolated_atom: 0.01
            liquid_high: 10.0
        contains:
            liquid: 5.0
    descriptors:
        - 2b:
            kernel_setups:
                - type: se
                  length_scale: 1.0
                  energy_scale: 10.0
                  n_sparse: 10.0
        ...
    run: false (true by default, data not specified => use all in context)
- branch:
    - modify_data:
        label: abc
        mod:
            - add_noise: {sigma, ...}
    - run_fit:
        fit_label: gap2b
        pot_label: gap2b_noise1
        data: abc
    - save:
        ... gap2b_noise1 ...
    - predict:
        (use current pot by default)
        on: abc
        to: abc_nv_1
        analyze: [rmse: {group_by: ..., save_to}, (elastic)]
 */

#endif
