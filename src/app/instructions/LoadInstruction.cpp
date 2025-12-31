#include "app/instructions/LoadInstruction.hpp"

#include <utility>

namespace jgap {

    std::shared_ptr<LoadInstruction> LoadInstruction::fromDataNode(const DataNode& params) {

        // Only the file name.
        if (params.type == DataNode::Type::STRING) {
            return std::make_shared<LoadInstruction>(params.asString());
        }

        if (params.contains("label")) {
            return std::make_shared<LoadInstruction>(
                REQUIRE(params, "file"),
                REQUIRE(params, "file").asString()
                );
        }

        return std::make_shared<LoadInstruction>(REQUIRE(params, "file"));
    }

    LoadInstruction::LoadInstruction(const std::string &file) {
        file_ = file;
        label_ = withoutExtension(file_);
    }

    LoadInstruction::LoadInstruction(const std::string& file, const std::string& label)
        : file_(file), label_(label) {
    }

    void LoadInstruction::execute(PipelineState &state) {
        if (file_.ends_with(".xyz")) {

            auto data = readXyz(file_);
            state.addData(label_, DataSetWithMeta{.data = data, .saved = true, .to_be_saved = false});

        } else if (file_.ends_with(".yaml")) {

            std::ifstream inFile(file_);
            if (!inFile.is_open()) {
                JGAP_LOG_ERROR("Could not open file {}", file_);
            }

            std::ostringstream instructionsStringStream;
            instructionsStringStream << inFile.rdbuf();

            YAML::Node yamlNode;
            try {
                yamlNode = YAML::Load(instructionsStringStream.str());
            } catch (std::exception& e) {
                JGAP_LOG_AND_THROW("Could not parse xml pot-file: {}", e.what());
            }

            DataNode potParams = YamlConverter::fromYaml(yamlNode);

            state.addPotential(label_, PotentialWithMeta{
                .potential = REGISTRY_GET(Potential, potParams),
                .saved = true
                });

        } else if (file_.ends_with(".xml")) {

            pugi::xml_document xml_doc;
            try {
                xml_doc.load_file(file_.c_str());
            } catch (std::exception& e) {
                JGAP_LOG_AND_THROW("Could not parse xml pot-file: {}", e.what());
            }

            state.addPotential(label_, PotentialWithMeta{
                .potential = QuipXmlConverter::transform(xml_doc.document_element()),
                .saved = false,
                .to_be_saved = true
                });
        }
    }
}