#ifndef JGAP_EXITINSTRUCTION_HPP
#define JGAP_EXITINSTRUCTION_HPP

#include "PipelineInstruction.hpp"

namespace jgap {
    class ExitInstruction : public PipelineInstruction {
    public:
        SETUP_PARSER(PipelineInstruction, ExitInstruction, exit)
        ExitInstruction();
        void execute(PipelineState &state) override;
    };
}

#endif
