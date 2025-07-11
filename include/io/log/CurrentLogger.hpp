#ifndef CURRENTLOGGER_HPP
#define CURRENTLOGGER_HPP

#include "Logger.hpp"
#include "StdoutLogger.hpp"

namespace jgap {
    class CurrentLogger {
    public:
        static shared_ptr<Logger> get() {
            return _logger;
        }
        static void setLogger(const shared_ptr<Logger> &logger) {
            _logger = logger;
        }
    private:
        inline static shared_ptr<Logger> _logger = make_shared<StdoutLogger>();
    };
}

#endif //CURRENTLOGGER_HPP
