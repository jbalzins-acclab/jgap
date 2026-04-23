#ifndef JGAP_BRANCHINSTRUCTION_HPP
#define JGAP_BRANCHINSTRUCTION_HPP

#include "PipelineInstruction.hpp"
#include "app/ProcessingPipeline.hpp"

namespace jgap {
    class BranchInstruction : public PipelineInstruction {
    public:
        SETUP_PARSER(PipelineInstruction, BranchInstruction, branch)

        inline static std::deque<size_t> LEVELS = {};

        BranchInstruction(std::unique_ptr<ProcessingPipeline> pipeline, const std::optional<std::string> &label = {})
            : pipeline_(std::move(pipeline)), label_(label) {};

        void execute(PipelineState &state) override;
    private:
        std::unique_ptr<ProcessingPipeline> pipeline_;
        std::optional<std::string> label_;
    };
}

#endif