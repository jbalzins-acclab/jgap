#include "app/instructions/LoadInstruction.hpp"

#include <utility>

namespace jgap {

    std::shared_ptr<LoadInstruction> LoadInstruction::fromDataNode(const DataNode& params) {

        // Only the file name.
        if (params.type == DataNode::Type::STRING) {
            return std::make_shared<LoadInstruction>(params.asString(), std::optional<std::string>{});
        }

        return std::make_shared<LoadInstruction>(
            require(params, "file"),
            params.getOrDefault("label", require(params, "file").asString())
            );
    }

    LoadInstruction::LoadInstruction(const std::string& file, const std::optional<std::string>& label)
        : file_(file)
    {
        if (label.has_value()) {
            label_ = label.value();
        } else {
            label_ = file_;
        }
    }

    void LoadInstruction::execute(PipelineState &state) {
        if (file_.ends_with(".xyz")) {

            assert(!state.data_sets.contains(label_) && "Data label already taken");
            state.data_sets[label_] = readXyz(file_);

        } else if (file_.ends_with(".yaml")) {

            assert(!state.potentials.contains(label_) && "Potential label already taken");

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

            state.potentials[label_] = PotentialWithMeta{
                .potential = REGISTRY_GET(Potential, potParams),
                .saved = true
                };

        } else if (file_.ends_with(".xml")) {

            assert(!state.potentials.contains(label_) && "Potential label already taken");

            pugi::xml_document xml_doc;
            try {
                xml_doc.load_file(file_.c_str());
            } catch (std::exception& e) {
                JGAP_LOG_AND_THROW("Could not parse xml pot-file: {}", e.what());
            }
            state.potentials[label_] = PotentialWithMeta{
                .potential = QuipXmlConverter::transform(xml_doc.document_element()),
                .saved = false
                };
        }
    }
}