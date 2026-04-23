#ifndef JGAP_LOADINSTRUCTION_HPP
#define JGAP_LOADINSTRUCTION_HPP

#include "PipelineInstruction.hpp"

#include <cassert>
#include <pugixml.hpp>

#include "io/convert/QuipXmlConverter.hpp"
#include "io/convert/YamlConverter.hpp"

namespace jgap {
    class LoadInstruction : public PipelineInstruction {
    public:
        SETUP_PARSER(PipelineInstruction, LoadInstruction, load)

        LoadInstruction(const std::string& file);
        LoadInstruction(const std::string& file, const std::string& label);
        void execute(PipelineState &state) override;

    private:
        std::string file_;
        std::string label_;
    };
}

#endif
