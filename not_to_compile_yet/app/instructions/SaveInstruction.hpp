#ifndef JGAP_SAVEINSTRUCTION_HPP
#define JGAP_SAVEINSTRUCTION_HPP

#include "PipelineInstruction.hpp"

namespace jgap {
    class SaveInstruction : public PipelineInstruction {
    public:
        SETUP_PARSER(PipelineInstruction, SaveInstruction, "save")
        SaveInstruction() : labels_({}) {}
        SaveInstruction(std::string label, std::string file) : labels_({label}) {}
        SaveInstruction(std::vector<std::string> labels);
        void execute(PipelineState &) override;

    private:
        std::string> labels_;
    };
}

#endif