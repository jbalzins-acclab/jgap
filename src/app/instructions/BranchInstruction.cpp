#include "BranchInstruction.hpp"

namespace jgap {
    std::shared_ptr<BranchInstruction> BranchInstruction::fromDataNode(const DataNode &params) {

        std::optional<std::string> label;

        if (params.type != DataNode::Type::OBJECT && params.contains("label")) {
            label = REQUIRE(params, "label").asString();
        } else if (params.type != DataNode::Type::ARRAY) {
            JGAP_LOG_AND_THROW("Branch instruction must be either an array or an object. Got: {}", params.toString());
        }

        auto branchPipeline = std::make_unique<ProcessingPipeline>(ProcessingPipeline::fromDataNode(params));

        return std::make_shared<BranchInstruction>(std::move(branchPipeline), label);
    }

    void BranchInstruction::execute(PipelineState &state) {
        pipeline_->changeState(state);

        try {
            pipeline_->execute();
        } catch (const std::exception& e) {
            if (label_.has_value()) {
                JGAP_LOG_ERROR("Branch {} failed: {}", label_.value(), e.what());
            } else {
                JGAP_LOG_ERROR("Branch failed: {}", e.what());
            }
        }
    }
}