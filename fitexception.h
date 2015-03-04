#ifndef ____EuiGhCCk
#define ____EuiGhCCk
#include <exception>
#include <string>
#include "fit_gen.h"
class FitException:public std::exception{
public:
	FitException(std::string msg);
	virtual ~FitException() throw(){}
	virtual const char* what() const throw();
private:
	std::string m_msg;
};
#endif
