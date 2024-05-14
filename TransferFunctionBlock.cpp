//
// Created by issa on 14/05/24.
//

#include "TransferFunctionBlock.h"
#include <iostream>
#include <sstream>

/**
 * @brief Pads a vector with zeros at the beginning to match a specified size.
 *
 * @param v The vector to be padded.
 * @param newSize The size to pad the vector to.
 */
void padWithZeros(Eigen::VectorXd &v, int newSize) {
    int originalSize = v.size();
    if (originalSize < newSize) {
        Eigen::VectorXd temp = Eigen::VectorXd::Zero(newSize);
        temp.segment(newSize - originalSize, originalSize) = v;
        v = temp;
    }
}

/**
 * @brief Adds two vectors element-wise, padding with zeros if necessary to match sizes.
 *
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return Eigen::VectorXd The result of adding the two vectors.
 */
Eigen::VectorXd addVectors(Eigen::VectorXd v1, Eigen::VectorXd v2) {
    int size1 = v1.size();
    int size2 = v2.size();
    int maxSize = std::max(size1, size2);

    // Pad the smaller vector with zeros at the beginning
    padWithZeros(v1, maxSize);
    padWithZeros(v2, maxSize);

    return v1 + v2;
}

/**
 * @brief Performs polynomial convolution of two vectors.
 *
 * @param v1 The first vector (representing polynomial coefficients).
 * @param v2 The second vector (representing polynomial coefficients).
 * @return Eigen::VectorXd The result of the convolution.
 */
Eigen::VectorXd convolve(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2) {
    int size1 = v1.size();
    int size2 = v2.size();
    int resultSize = size1 + size2 - 1;

    Eigen::VectorXd result = Eigen::VectorXd::Zero(resultSize);

    for (int i = 0; i < size1; ++i) {
        for (int j = 0; j < size2; ++j) {
            result[i + j] += v1[i] * v2[j];
        }
    }

    return result;
}

/**
 * @brief Formats a polynomial from its coefficients into a string.
 *
 * @param coefficients The coefficients of the polynomial.
 * @return std::string The formatted polynomial.
 */
std::string formatPolynomial(const Eigen::VectorXd &coefficients) {
    std::ostringstream oss;
    bool first = true;

    for (int i = 0; i < coefficients.size(); ++i) {
        double coeff = coefficients[i];
        int power = coefficients.size() - 1 - i;

        if (coeff != 0) {
            if (!first && coeff > 0) {
                oss << " + ";
            } else if (coeff < 0) {
                oss << " - ";
                coeff = -coeff;  // Make coefficient positive for printing
            }

            if (coeff != 1 || power == 0) {
                oss << coeff;
            }

            if (power > 0) {
                oss << "s";
                if (power > 1) {
                    oss << "^" << power;
                }
            }

            first = false;
        }
    }

    if (first) {
        // If all coefficients are zero, the polynomial is 0
        oss << "0";
    }

    return oss.str();
}

/**
 * @brief Default constructor initializing a transfer function block with numerator 0 and denominator 1.
 */
TransferFunctionBlock::TransferFunctionBlock()
        : numerator_(Eigen::VectorXd::Zero(1)), denominator_(Eigen::VectorXd::Ones(1)) {}

/**
 * @brief Constructor initializing a transfer function block with given numerator and denominator.
 *
 * @param numerator The numerator coefficients of the transfer function.
 * @param denominator The denominator coefficients of the transfer function.
 */
TransferFunctionBlock::TransferFunctionBlock(const Eigen::VectorXd &numerator, const Eigen::VectorXd &denominator)
        : numerator_(numerator), denominator_(denominator) {
}

/**
 * @brief Prints the transfer function in a readable format.
 */
void TransferFunctionBlock::PrintTransferFunction() const {
    std::string num_str = formatPolynomial(numerator_);
    std::string den_str = formatPolynomial(denominator_);

    std::cout << "Transfer Function: " << std::endl;
    std::cout << "      " << num_str << std::endl;
    std::cout << "H(s) = -----------------" << std::endl;
    std::cout << "      " << den_str << std::endl;
}

/**
 * @brief Connects two transfer function blocks in series.
 *
 * @param block1 The first transfer function block.
 * @param block2 The second transfer function block.
 * @return TransferFunctionBlock The resulting transfer function block from the series connection.
 */
TransferFunctionBlock TransferFunctionBlock::SeriesConnection(const TransferFunctionBlock &block1, const TransferFunctionBlock &block2) {
    Eigen::VectorXd new_num = convolve(block1.numerator_, block2.numerator_);
    Eigen::VectorXd new_den = convolve(block1.denominator_, block2.denominator_);
    return TransferFunctionBlock{new_num, new_den};
}

/**
 * @brief Connects two transfer function blocks in parallel.
 *
 * @param block1 The first transfer function block.
 * @param block2 The second transfer function block.
 * @return TransferFunctionBlock The resulting transfer function block from the parallel connection.
 */
TransferFunctionBlock TransferFunctionBlock::ParallelConnection(const TransferFunctionBlock &block1, const TransferFunctionBlock &block2) {
    // Compute the resulting numerator and denominator polynomials through addition
    Eigen::VectorXd new_num = addVectors(block1.numerator_, block2.numerator_);
    Eigen::VectorXd new_den = block1.denominator_;

    // Return a new TransferFunctionBlock with the resulting coefficients
    return TransferFunctionBlock{new_num, new_den};
}

/**
 * @brief Applies feedback to a transfer function block.
 *
 * @param block1 The transfer function block representing the forward path.
 * @param block2 The transfer function block representing the feedback path.
 * @param feedback_gain The feedback gain.
 * @return TransferFunctionBlock The resulting transfer function block with feedback applied.
 */
TransferFunctionBlock TransferFunctionBlock::FeedbackConnection(const TransferFunctionBlock &block1, const TransferFunctionBlock &block2, double feedback_gain) {
    // Convolve block1 numerator with block2 denominator
    Eigen::VectorXd new_num = convolve(block1.numerator_, block2.denominator_);
    // Convolve block1 denominator with block2 denominator
    Eigen::VectorXd block1_den_block2_den_convolved = convolve(block1.denominator_, block2.denominator_);
    // Convolve block1 numerator with block2 numerator
    Eigen::VectorXd block1_num_block2_num_convolved = convolve(block1.numerator_, block2.numerator_);
    // Scale the convolution result by the feedback gain
    block1_num_block2_num_convolved *= feedback_gain;
    // The new denominator
    Eigen::VectorXd new_den = addVectors(block1_num_block2_num_convolved, block1_den_block2_den_convolved);
    // Return the resulting TransferFunctionBlock
    return TransferFunctionBlock{new_num, new_den};
}
