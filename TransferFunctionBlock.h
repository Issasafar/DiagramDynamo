//
// Created by issa on 14/05/24.
//

#ifndef DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H
#define DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H

#include <vector>
#include <deque>
#include <Eigen/Dense>

/**
 * @class TransferFunctionBlock
 * @brief This class represents a transfer function block in a control system.
 *
 * The TransferFunctionBlock class provides methods for creating and manipulating transfer functions.
 * It supports operations such as printing the transfer function, and connecting transfer functions in series,
 * parallel, and feedback configurations.
 */
class TransferFunctionBlock {
public:
    /**
     * @brief Default constructor for TransferFunctionBlock.
     */
    TransferFunctionBlock();

    /**
     * @brief Parameterized constructor for TransferFunctionBlock.
     * @param numerator The numerator coefficients of the transfer function.
     * @param denominator The denominator coefficients of the transfer function.
     */
    TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator);

    /**
     * @brief Prints the transfer function in a readable format.
     */
    void PrintTransferFunction() const;

    /**
     * @brief Creates a new TransferFunctionBlock by connecting two blocks in series.
     * @param block1 The first TransferFunctionBlock.
     * @param block2 The second TransferFunctionBlock.
     * @return A new TransferFunctionBlock representing the series connection of block1 and block2.
     */
    static TransferFunctionBlock SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);

    /**
     * @brief Creates a new TransferFunctionBlock by connecting two blocks in parallel.
     * @param block1 The first TransferFunctionBlock.
     * @param block2 The second TransferFunctionBlock.
     * @return A new TransferFunctionBlock representing the parallel connection of block1 and block2.
     */
    static TransferFunctionBlock ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);

    /**
     * @brief Creates a new TransferFunctionBlock representing the application of feedback to an existing block.
     * @param block The TransferFunctionBlock object to which feedback is applied.
     * @param feedback_block The TransferFunctionBlock object providing the feedback.
     * @param feedback_gain The feedback gain value.
     * @return A new TransferFunctionBlock object representing the system with feedback applied.
     */
    static TransferFunctionBlock FeedbackConnection(const TransferFunctionBlock& block, const TransferFunctionBlock& feedback_block, double feedback_gain);

private:
    Eigen::VectorXd numerator_; ///< The numerator coefficients of the transfer function.
    Eigen::VectorXd denominator_; ///< The denominator coefficients of the transfer function.
};

#endif //DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H
