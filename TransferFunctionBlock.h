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
    // Enum for signal types (currently unused, can be added for future expansion)
    // enum class SignalType { CONTINUOUS, DISCRETE };

    /**
     * Default constructor that initializes both numerator and denominator with all coefficients set to zero (representing a constant zero system).
     */
    TransferFunctionBlock();

    /**
     * Constructor that takes separate numerator and denominator coefficients (Eigen vectors) for defining custom transfer functions.
     * @param numerator Eigen vector representing the numerator coefficients.
     * @param denominator Eigen vector representing the denominator coefficients.
     */
    TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator);

    /**
     * Sends a new input to the system and returns the calculated output (continuous time).
     * @param input The new input value to be fed into the system.
     * @return The calculated output of the system based on the given input and past history.
     */
    double SendInput(double input);

    /**
     * Placeholder function intended to provide access to the current output of the system based on past inputs/outputs.
     * Needs to be implemented based on specific requirements (e.g., retrieving the most recent output for systems with delays).
     * @return (Currently a placeholder, modify for your specific needs)
     */
    double GetOutput() const;

    // Static methods for series, parallel, and feedback connections

    /**
     * Creates a new TransferFunctionBlock object representing the series connection of two existing blocks.
     * @param block1 The first TransferFunctionBlock object.
     * @param block2 The second TransferFunctionBlock object.
     * @return A new TransferFunctionBlock object representing the series connection of the two input blocks.
     */
    static TransferFunctionBlock SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);

    /**
     * Creates a new TransferFunctionBlock object representing the parallel connection of two existing blocks.
     * @param block1 The first TransferFunctionBlock object.
     * @param block2 The second TransferFunctionBlock object.
     * @return A new TransferFunctionBlock object representing the parallel connection of the two input blocks.
     */
    static TransferFunctionBlock ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2);

    /**
     * Creates a new TransferFunctionBlock object representing the application of feedback to an existing block.
     * @param block The TransferFunctionBlock object to which feedback is applied.
     * @param feedback_gain The feedback gain value.
     * @return A new TransferFunctionBlock object representing the system with feedback applied.
     */
    static TransferFunctionBlock FeedbackConnection(const TransferFunctionBlock& block, double feedback_gain);

private:
    // Internal state variables
    Eigen::VectorXd numerator_;
    Eigen::VectorXd denominator_;
    std::deque<double> pastInputs_;
    std::deque<double> pastOutputs_;

    // Function to evaluate the polynomial (using Eigen)
    double EvaluatePolynomial(const Eigen::VectorXd& p, double x) const;
};

#endif //DIAGRAMDYNAMO_TRANSFERFUNCTIONBLOCK_H
