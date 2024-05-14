//
// Created by issa on 14/05/24.
//

#include "TransferFunctionBlock.h"
#include "TransferFunctionBlock.h"

TransferFunctionBlock::TransferFunctionBlock() : numerator_(Eigen::VectorXd::Zero(1)), denominator_(Eigen::VectorXd::Ones(1)) {}

TransferFunctionBlock::TransferFunctionBlock(const Eigen::VectorXd& numerator, const Eigen::VectorXd& denominator) :
        numerator_(numerator), denominator_(denominator) {}

double TransferFunctionBlock::SendInput(double input) {
    pastInputs_.push_front(input);
    if (pastInputs_.size() > numerator_.size()) {
        pastInputs_.pop_back();
    }

    double output = 0.0;
    // Calculate output using Eigen for polynomial evaluation (numerator and denominator)
    for (int i = 0; i < numerator_.size(); ++i) {
        output += EvaluatePolynomial(numerator_, pastInputs_[i]);
    }
    for (int i = 0; i < denominator_.size(); ++i) {
        if (i > 0) {
            output -= EvaluatePolynomial(denominator_, pastOutputs_[i-1]);
        }
    }

    pastOutputs_.push_front(output);
    if (pastOutputs_.size() > denominator_.size()) {
        pastOutputs_.pop_back();
    }

    return output;
}

double TransferFunctionBlock::GetOutput() const {
    // Implement logic to return current output based on past inputs/outputs (if applicable)
    return 0.0; // Placeholder, modify for your specific needs
}


double EvaluatePolynomial(const Eigen::VectorXd& coefficients, double x) const {
    double result = 0.0;
    double x_power = 1.0;  // Start with x^0 which is 1

    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * x_power;
        x_power *= x;  // Update x_power to x^(i+1)
    }

    return result;
}

// Static methods for connection

TransferFunctionBlock TransferFunctionBlock::SeriesConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    Eigen::VectorXd new_num = block1.numerator_ * block2.numerator_;
    Eigen::VectorXd new_den = block1.denominator_ * block2.denominator_;
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::ParallelConnection(const TransferFunctionBlock& block1, const TransferFunctionBlock& block2) {
    Eigen::VectorXd new_num = block1.numerator_ * block2.denominator_ +
                              block2.numerator_ * block1.denominator_;
    Eigen::VectorXd new_den = block1.denominator_ * block2.denominator_;
    return TransferFunctionBlock(new_num, new_den);
}

TransferFunctionBlock TransferFunctionBlock::FeedbackConnection(const TransferFunctionBlock& block, double feedback_gain) {
    Eigen::VectorXd new_num = block.numerator_ * (1.0 + feedback_gain);
    Eigen::VectorXd new_den = block.denominator_ + feedback_gain * block.denominator_;
    return TransferFunctionBlock(new_num, new_den);
}

