# filtration

These documents with .cpp and .py are used for simulating the filtration process.

testparticlefiltration.cpp could be used for the normal slurry with small particles.

testflocfiltration.cpp could be used for slurry after flocculation with big floc.

filter.cpp can be used for testing the initial filter.

All the .py documents are used for analysis the simulation.

Before using these .cpp documents, you need to download the documents named "dlvo" forked from hhucyl/dlvo,
and copy these .cpp into the dlvo.

These simulation code has not been finished yet, I still work on it and will update these codes recently.

Any questions can be sent to the Email address: wusilinhhu@126.com  or 1070195170@qq.com.

Be happy with simulation.

Cheers.

2019-07-18-14:05

updata record:

{

Three documents have been updata:

testfilter1.0.cpp                  //this is for test the filter;

filter1.0.py                       //this is for analysis the filter;

10.Filter渗透系数的误差测试报告       // this is some analysis for this code.

}

2019-07-22-16:10

updata record:

{

Two documents have been updata:

testflocfiltration1.01.cpp         //updata the supplementary floc code

filtrationtest1.0.py              //this is for analysising the test of filtration

}

2019-07-23-18:09

updata record:

{

Two documents have been updata:

testflocfiltration1.02.cpp         //updata the Gn and MU

testparticlefiltration1.0         //updata the Gn and MU

}

2019-08-02-11:32

updata record:

{

the newest one for floc filtration test is testflocfiltration1.02.cpp;

the newest one for particle filtration test is testparticlefiltration1.0.cpp;

the newest one for filter test is testfilter1.0.cpp;

the newest one for particle or floc filtration test analysis is filtration1.04.py;

the newest one for filter test analysis is filter1.0.py;

}

2019-08-04-12:13

updata record:

{

after huge testing on testflocfiltration1.02.cpp, we found there was small bug on the floc supplement, some floc will be located beyond the ny area, that lead some bug, the floc supplement will be slower or broken when entering into the nx*ny area. so we fix this bug, and update the testflocfiltration1.03.cpp 

}

2019-08-04-12:13

updata record:

{

the newest one for floc filtration test is testflocfiltration1.03.cpp;

the newest one for particle filtration test is testparticlefiltration1.0.cpp;

the newest one for filter test is testfilter1.0.cpp;

the newest one for particle or floc filtration test analysis is filtration1.04.py;

the newest one for filter test analysis is filter1.0.py;

}

2019-08-08-11:51

updata record:

{

the newest one for floc filtration test is testflocfiltration1.03.cpp;

the newest one for particle filtration test is testparticlefiltration1.01.cpp;

the newest one for filter test is testfilter1.0.cpp;

the newest one for particle or floc filtration test analysis is filtration1.04.py;

the newest one for filter test analysis is filter1.0.py;

}
