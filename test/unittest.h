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
