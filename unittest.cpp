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
		printf("\tTest : %s : ...", (*iter)->GetDescription().c_str());
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
