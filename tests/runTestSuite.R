
testsuite.bootES = RUnit::defineTestSuite(
                             "bootES",
                             dirs=system.file("tests", package="bootES"),
                             testFileRegexp="^runit",
                             testFuncRegexp="^test")

testResult = RUnit::runTestSuite(testsuite.bootES)
