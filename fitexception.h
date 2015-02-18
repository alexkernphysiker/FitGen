#ifndef APROEXCEPTION_H
#define APROEXCEPTION_H
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
#endif // APROEXCEPTION_H
