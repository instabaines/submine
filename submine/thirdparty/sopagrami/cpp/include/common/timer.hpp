#pragma once
#include <chrono>

namespace common {

class Timer {
public:
    using clock = std::chrono::steady_clock;
    void start() { t0_ = clock::now(); }
    double ms() const {
        auto t1 = clock::now();
        return std::chrono::duration<double, std::milli>(t1 - t0_).count();
    }
private:
    clock::time_point t0_;
};

} // namespace common
