#include <iostream>
#include "TransferFunctionBlock.h"

int main() {
    // Define numerator and denominator coefficients for the transfer function
    Eigen::VectorXd numerator(3);
    numerator << 1.0, 0.5, 0.2; // Example coefficients for the numerator (1 + 0.5s + 0.2s^2)

    Eigen::VectorXd denominator(3);
    denominator << 1.0, -0.3, 0.1; // Example coefficients for the denominator (1 - 0.3s + 0.1s^2)

    // Create a TransferFunctionBlock object
    TransferFunctionBlock tfBlock(numerator, denominator);

    // Simulate sending inputs to the transfer function block
    std::vector<double> inputs = {1.0, 2.0, 3.0, 4.0, 5.0}; // Example input sequence
    for (double input : inputs) {
        double output = tfBlock.SendInput(input);
        std::cout << "Input: " << input << " -> Output: " << output << std::endl;
    }

    // Example of getting the current output (optional)
    double currentOutput = tfBlock.GetOutput();
    std::cout << "Current Output: " << currentOutput << std::endl;

    // Example of creating a series connection of two transfer function blocks
    TransferFunctionBlock block1(numerator, denominator);
    TransferFunctionBlock block2(numerator, denominator);
    TransferFunctionBlock seriesBlock = TransferFunctionBlock::SeriesConnection(block1, block2);

    // Example of creating a parallel connection of two transfer function blocks
    TransferFunctionBlock parallelBlock = TransferFunctionBlock::ParallelConnection(block1, block2);

    // Example of creating a feedback connection of a transfer function block
    double feedback_gain = 0.5;
    TransferFunctionBlock feedbackBlock = TransferFunctionBlock::FeedbackConnection(block1, feedback_gain);

    // Simulate sending inputs to the seriesBlock, parallelBlock, and feedbackBlock (optional)
    for (double input : inputs) {
        double seriesOutput = seriesBlock.SendInput(input);
        double parallelOutput = parallelBlock.SendInput(input);
        double feedbackOutput = feedbackBlock.SendInput(input);
        std::cout << "Series Connection Output: " << seriesOutput << std::endl;
        std::cout << "Parallel Connection Output: " << parallelOutput << std::endl;
        std::cout << "Feedback Connection Output: " << feedbackOutput << std::endl;
    }

    return 0;
}
