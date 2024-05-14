#include "TransferFunctionBlock.h"
#include <iostream>

int main() {
    // Define some example transfer function blocks
    Eigen::VectorXd num1(3);
    num1 << 1, 2, 3;
    Eigen::VectorXd den1(2);
    den1 << 1, 1;
    TransferFunctionBlock block1(num1, den1);

    Eigen::VectorXd num2(2);
    num2 << 4, 5;
    Eigen::VectorXd den2(2);
    den2 << 1, 2;
    TransferFunctionBlock block2(num2, den2);

    Eigen::VectorXd num3(3);
    num3 << 1, 0, 1;
    Eigen::VectorXd den3(3);
    den3 << 1, 3, 2;
    TransferFunctionBlock block3(num3, den3);

    Eigen::VectorXd num4(2);
    num4 << 2, 3;
    Eigen::VectorXd den4(3);
    den4 << 1, 4, 4;
    TransferFunctionBlock block4(num4, den4);

    Eigen::VectorXd num5(1);
    num5 << 1;
    Eigen::VectorXd den5(2);
    den5 << 1, 2;
    TransferFunctionBlock block5(num5, den5);

    // Print the initial blocks
    std::cout << "Block 1:" << std::endl;
    block1.PrintTransferFunction();
    std::cout << std::endl;

    std::cout << "Block 2:" << std::endl;
    block2.PrintTransferFunction();
    std::cout << std::endl;

    std::cout << "Block 3:" << std::endl;
    block3.PrintTransferFunction();
    std::cout << std::endl;

    std::cout << "Block 4:" << std::endl;
    block4.PrintTransferFunction();
    std::cout << std::endl;

    std::cout << "Block 5:" << std::endl;
    block5.PrintTransferFunction();
    std::cout << std::endl;

    // Example 1: Series connection of block1 and block2
    TransferFunctionBlock series1 = TransferFunctionBlock::SeriesConnection(block1, block2);
    std::cout << "Series Connection of block1 and block2:" << std::endl;
    series1.PrintTransferFunction();
    std::cout << std::endl;

    // Example 2: Parallel connection of block3 and block4
    TransferFunctionBlock parallel1 = TransferFunctionBlock::ParallelConnection(block3, block4);
    std::cout << "Parallel Connection of block3 and block4:" << std::endl;
    parallel1.PrintTransferFunction();
    std::cout << std::endl;

    // Example 3: Feedback connection of block1 and block3 with gain 1
    TransferFunctionBlock feedback1 = TransferFunctionBlock::FeedbackConnection(block1, block3, 1.0);
    std::cout << "Feedback Connection of block1 and block3 with gain 1:" << std::endl;
    feedback1.PrintTransferFunction();
    std::cout << std::endl;

    // Example 4: Series connection of (block1 and block2) and block4
    TransferFunctionBlock series2 = TransferFunctionBlock::SeriesConnection(series1, block4);
    std::cout << "Series Connection of (block1 and block2) and block4:" << std::endl;
    series2.PrintTransferFunction();
    std::cout << std::endl;

    // Example 5: Feedback connection of (block1 in series with block2) and (block3 in parallel with block4) with gain 2
    TransferFunctionBlock series3 = TransferFunctionBlock::SeriesConnection(block1, block2);
    TransferFunctionBlock parallel2 = TransferFunctionBlock::ParallelConnection(block3, block4);
    TransferFunctionBlock feedback2 = TransferFunctionBlock::FeedbackConnection(series3, parallel2, 2.0);
    std::cout << "Feedback Connection of (block1 in series with block2) and (block3 in parallel with block4) with gain 2:" << std::endl;
    feedback2.PrintTransferFunction();
    std::cout << std::endl;

    return 0;
}
