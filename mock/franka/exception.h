// Mock of <franka/exception.h> (libfranka 0.9.2). See mock/franka/README.md.
#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "franka/log.h"

namespace franka {

    struct Exception : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    struct ModelException : public Exception {
        using Exception::Exception;
    };

    struct NetworkException : public Exception {
        using Exception::Exception;
    };

    struct ProtocolException : public Exception {
        using Exception::Exception;
    };

    struct IncompatibleVersionException : public Exception {
        IncompatibleVersionException(uint16_t server, uint16_t library) noexcept
            : Exception("libfranka: Incompatible library version"), server_version(server), library_version(library)
        {
        }
        const uint16_t server_version;
        const uint16_t library_version;
    };

    struct ControlException : public Exception {
        explicit ControlException(const std::string& what, std::vector<franka::Record> records = {}) noexcept
            : Exception(what), log(std::move(records))
        {
        }
        const std::vector<franka::Record> log;
    };

    struct CommandException : public Exception {
        using Exception::Exception;
    };

    struct RealtimeException : public Exception {
        using Exception::Exception;
    };

    struct InvalidOperationException : public Exception {
        using Exception::Exception;
    };

}  // namespace franka
