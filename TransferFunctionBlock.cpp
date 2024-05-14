#include "TransferFunctionBlock.h"
#include <iostream>
#include <sstream>

// Helper function to perform polynomial convolution
Eigen::VectorXd convolve(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(a.size() + b.size() - 1);
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < b.size(); ++j) {
            result[i + j] += a[i] * b[j];
        }
    }
    return result;
}

// Helper function to format a polynomial from coefficients
std::string formatPolynomial(const Eigen::VectorXd& coefficients) {
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

TransferFunctionBlock::TransferFunctionBlock()
        : numerator_(Eigen::VectorXd::Zero(1)), denominator_(Eigen::VectorXd::Ones(1)) {}

TransferFunctionBlock::TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator)
        : numerator_(numerator), denominator_(denominator) {
}



void TransferFunctionBlock::PrintTransferFunction() const {
    std::string num_str = formatPolynomial(numerator_);
    std::string den_str = formatPolynomial(denominator_);

    std::cout << "Transfer Function: " << std::endl;
    std::cout << "      " << num_str << std::endl;
    std::cout << "H(s) = -----------------" << std::endl;
    std::cout << "      " << den_str << std::endl;
}

TransferFunctionBlock TransferFunctionBlock::SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    Eigen::VectorXd new_num = convolve(block1.numerator_, block2.numerator_);
    Eigen::VectorXd new_den = convolve(block1.denominator_, block2.denominator_);
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    // Compute the resulting numerator and denominator polynomials through addition
    Eigen::VectorXd new_num = block1.numerator_ + block2.numerator_;
    Eigen::VectorXd new_den = block1.denominator_;

    // Return a new TransferFunctionBlock with the resulting coefficients
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::FeedbackConnection(const TransferFunctionBlock& gs, const TransferFunctionBlock& fs, double k) {
    std::cout << "Sizes: gs numerator: " << gs.numerator_.size() << ", gs denominator: " << gs.denominator_.size()
              << ", fs numerator: " << fs.numerator_.size() << ", fs denominator: " << fs.denominator_.size() << std::endl;

    // Convolve gs numerator with fs denominator
    Eigen::VectorXd gs_num_fs_den_convolved = convolve(gs.numerator_, fs.denominator_);
    std::cout << "Size of gs_num_fs_den_convolved: " << gs_num_fs_den_convolved.size() << std::endl;

    // Scale the convolution by the feedback gain
    gs_num_fs_den_convolved *= k;

    // Convolve gs denominator with fs numerator
    Eigen::VectorXd gs_den_fs_num_convolved = convolve(gs.denominator_, fs.numerator_);
    std::cout << "Size of gs_den_fs_num_convolved: " << gs_den_fs_num_convolved.size() << std::endl;

    // Convolve gs denominator with fs denominator
    Eigen::VectorXd gs_den_fs_den_convolved = convolve(gs.denominator_, fs.denominator_);
    std::cout << "Size of gs_den_fs_den_convolved: " << gs_den_fs_den_convolved.size() << std::endl;

    // Create the new numerator
    Eigen::VectorXd new_num = gs_num_fs_den_convolved;

    // Create the new denominator
    Eigen::VectorXd new_den = gs_den_fs_den_convolved + gs_den_fs_num_convolved;

    // Return the resulting TransferFunctionBlock
    return TransferFunctionBlock(new_num, new_den);
}

//TransferFunctionBlock TransferFunctionBlock::FeedbackConnection(const TransferFunctionBlock& gs, const TransferFunctionBlock& fs, double k) {
//    // Convolve gs denominator with fs numerator
//    Eigen::VectorXd gs_den_fs_num_convolved = convolve(gs.denominator_, fs.numerator_);
//
//    // Scale the convolution by the feedback gain
//    gs_den_fs_num_convolved *= k;
//
//    // Convolve gs numerator with fs denominator
//    Eigen::VectorXd gs_num_fs_den_convolved = convolve(gs.numerator_, fs.denominator_);
//
//    // Convolve gs denominator with fs denominator
//    Eigen::VectorXd gs_den_fs_den_convolved = convolve(gs.denominator_, fs.denominator_);
//
//    // Create the new numerator
//    Eigen::VectorXd new_num = gs_num_fs_den_convolved;
//
//    // Create the new denominator
//    Eigen::VectorXd new_den = gs_den_fs_den_convolved + gs_den_fs_num_convolved;
//
//    // Return the resulting TransferFunctionBlock
//    return TransferFunctionBlock(new_num, new_den);
//}


