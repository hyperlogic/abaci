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

#include "unittest.h"
#include <stdio.h>

TestSuite::~TestSuite()
{
	for(TestCaseVec::iterator iter = m_testCaseVec.begin(); iter != m_testCaseVec.end(); ++iter)
		delete *iter;
}
	
void TestSuite::AddTest(TestCase* testCase)
{
	m_testCaseVec.push_back(testCase);
}

bool TestSuite::RunTests() const
{
	printf("Running Test Suite : %s\n", m_description.c_str());
	bool suiteResult = true;
	for(TestCaseVec::const_iterator iter = m_testCaseVec.begin(); iter != m_testCaseVec.end(); ++iter)
	{
		printf("\tTest : %s : ... ", (*iter)->GetDescription().c_str());
		bool result = (*iter)->Test();
		suiteResult &= result;

		if (result)
			printf(" Success\n");
		else
			printf(" Failure\n");
	}
	
	printf("Test Suite %s Result : %s\n", m_description.c_str(), suiteResult ? "Success" : "Failure");
	return suiteResult;
}
