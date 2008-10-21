/*
The MIT License

Copyright (c) 2008 Anthony J. Thibault

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef UNITTEST_H
#define UNITTEST_H

#include <vector>
#include <string>

class TestCase
{
public:
	TestCase(const std::string& description) { m_description = description; }
	virtual ~TestCase() {};
	virtual bool Test() const = 0;

	const std::string& GetDescription() const { return m_description; }
protected:
	std::string m_description;
};

class TestSuite
{
public:
	TestSuite(const std::string& description) { m_description = description; }
	~TestSuite();
	
	void AddTest(TestCase* testCase);
	bool RunTests() const;
protected:
	std::string m_description;
	typedef std::vector<TestCase*> TestCaseVec;
	TestCaseVec m_testCaseVec;
};

#endif
