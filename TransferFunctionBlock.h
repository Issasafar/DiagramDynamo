//
// Created by issa on 14/05/24.
//

#ifndef DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H
#define DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H

#include <vector>
#include <deque>
#include <Eigen/Dense>

class TransferFunctionBlock {
public:
    TransferFunctionBlock();
    TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator);

    void PrintTransferFunction() const;

    static TransferFunctionBlock SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);
    static TransferFunctionBlock ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);
    /**
       * Creates a new TransferFunctionBlock object representing the application of feedback to an existing block.
       * @param block The TransferFunctionBlock object to which feedback is applied.
       * @param feedback_block The TransferFunctionBlock object providing the feedback.
       * @param feedback_gain The feedback gain value.
       * @return A new TransferFunctionBlock object representing the system with feedback applied.
       */
    static TransferFunctionBlock FeedbackConnection(const TransferFunctionBlock& block, const TransferFunctionBlock& feedback_block, double feedback_gain);

private:
    Eigen::VectorXd numerator_;
    Eigen::VectorXd denominator_;
};

#endif //DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H
