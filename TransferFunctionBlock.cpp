#include "TransferFunctionBlock.h"
#include <cassert>

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

TransferFunctionBlock::TransferFunctionBlock()
        : numerator_(Eigen::VectorXd::Zero(1)), denominator_(Eigen::VectorXd::Ones(1)) {}

TransferFunctionBlock::TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator)
        : numerator_(numerator), denominator_(denominator) {
    // Initialize past inputs and outputs
    pastInputs_ = std::deque<double>(numerator.size(), 0.0);
    pastOutputs_ = std::deque<double>(denominator.size() - 1, 0.0);
}

double TransferFunctionBlock::SendInput(double input) {
    pastInputs_.push_front(input);
    if (pastInputs_.size() > numerator_.size()) {
        pastInputs_.pop_back();
    }

    double output = 0.0;

    // Calculate the numerator part
    for (int i = 0; i < numerator_.size(); ++i) {
        if (i < pastInputs_.size()) {
            output += numerator_[i] * pastInputs_[i];
        }
    }

    // Calculate the denominator part (skip the first element which is the current output)
    for (int i = 1; i < denominator_.size(); ++i) {
        if (i - 1 < pastOutputs_.size()) {
            output -= denominator_[i] * pastOutputs_[i - 1];
        }
    }

    pastOutputs_.push_front(output);
    if (pastOutputs_.size() > denominator_.size() - 1) {
        pastOutputs_.pop_back();
    }

    return output;
}

double TransferFunctionBlock::GetOutput() const {
    // Return the most recent output
    if (!pastOutputs_.empty()) {
        return pastOutputs_.front();
    }
    return 0.0; // Return a default value if no output is available
}

double TransferFunctionBlock::EvaluatePolynomial(const Eigen::VectorXd& coefficients, double x) const {
    double result = 0.0;
    double x_power = 1.0; // Start with x^0 which is 1

    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * x_power;
        x_power *= x; // Update x_power to x^(i+1)
    }

    return result;
}

TransferFunctionBlock TransferFunctionBlock::SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    Eigen::VectorXd new_num = convolve(block1.numerator_, block2.numerator_);
    Eigen::VectorXd new_den = convolve(block1.denominator_, block2.denominator_);
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    Eigen::VectorXd new_num = convolve(block1.numerator_, block2.denominator_) + convolve(block2.numerator_, block1.denominator_);
    Eigen::VectorXd new_den = convolve(block1.denominator_, block2.denominator_);
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::FeedbackConnection(const TransferFunctionBlock& block, double feedback_gain) {
    Eigen::VectorXd new_num = block.numerator_;
    Eigen::VectorXd new_den = block.denominator_ + feedback_gain * convolve(block.numerator_, Eigen::VectorXd::Ones(1));
    return TransferFunctionBlock(new_num, new_den);
}
