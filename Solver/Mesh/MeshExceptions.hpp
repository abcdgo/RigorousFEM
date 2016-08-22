#ifndef _MESHEXCEPTIONS_H_
#define _MESHEXCEPTIONS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

template <typename ScalarType>
class BoundaryException : public std::runtime_error
{
    public:
    std::string message;

    BoundaryException(const std::string &info, const ScalarType &left, const ScalarType &right) : std::runtime_error(info)
    {
        std::ostringstream d;
        d << std::runtime_error::what() << "\n";
        d << "Left bound: " << left << "\n";
        d << "Right bound: " << right << "\n";
        message = d.str();
    }

    ~BoundaryException() throw() {}

    const char* what() const throw()
    {
        return message.c_str();
    }
};
#endif // MESHEXCEPTIONS_H_
