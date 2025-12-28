#ifndef JGAP_PIPELINEINSTRUCTION_HPP
#define JGAP_PIPELINEINSTRUCTION_HPP

#include "app/PipelineState.hpp"

namespace jgap {
    class PipelineInstruction {
    public:
        virtual ~PipelineInstruction() = default;
        virtual void execute(PipelineState&) = 0;
    };
}

#endif
