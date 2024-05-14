#include <iostream>
#include "TransferFunctionBlock.h"

using coef = Eigen::VectorXd;

int main() {
    // Define numerator and denominator coefficients for the transfer functions
    Eigen::VectorXd numerator(2);
    numerator << 1.0, 1.0; // Example coefficients for the numerator (s + 1)

    Eigen::VectorXd denominator(1);
    denominator << 1.0; // Example coefficients for the denominator (1)

    coef num(2);
    num << 1.0, 0.0; // Example coefficients for the numerator (s)

    coef den(1);
    den << 1.0; // Example coefficients for the denominator (1)

    // Create TransferFunctionBlock objects
    TransferFunctionBlock tfBlock(numerator, denominator);
    TransferFunctionBlock tf2Block(num, den);

    // Print the transfer functions
    std::cout << "Transfer Function 1: " << std::endl;
    tfBlock.PrintTransferFunction();
    std::cout << std::endl;

    std::cout << "Transfer Function 2: " << std::endl;
    tf2Block.PrintTransferFunction();
    std::cout << std::endl;

    // Perform feedback connection and print the resulting transfer function
    TransferFunctionBlock feedbackResult = TransferFunctionBlock::FeedbackConnection(tfBlock, tf2Block, -1);
    std::cout << "Resulting Transfer Function after Feedback Connection: " << std::endl;
    feedbackResult.PrintTransferFunction();

    return 0;
}
