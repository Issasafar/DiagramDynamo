#include <iostream>
#include "TransferFunctionBlock.h"

int main() {
    // Define numerator and denominator coefficients for the transfer function
    Eigen::VectorXd numerator(3);
    numerator << 0.2, 0.5, 1.0; // Example coefficients for the numerator (0.2 + 0.5s + s^2)

    Eigen::VectorXd denominator(3);
    denominator << 1.0, -0.3, 0.1; // Example coefficients for the denominator (1 - 0.3s + 0.1s^2)

    // Create a TransferFunctionBlock object
    TransferFunctionBlock tfBlock(numerator, denominator);

    // Print the transfer function
    tfBlock.PrintTransferFunction();

    // Simulate sending inputs to the transfer function block
    std::vector<double> inputs = {1.0, 2.0, 3.0, 4.0, 5.0}; // Example input sequence
    for (double input : inputs) {
        double output = tfBlock.SendInput(input);
        std::cout << "Input: " << input << " -> Output: " << output << std::endl;
    }

    return 0;
}
