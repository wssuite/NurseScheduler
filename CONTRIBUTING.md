# CONTRIBUTING

Any contribution needs to be merged in the master with a pull request.

Once the pull request is created, github action will run two workflows:
 - code-analysis: cpplint is run on the code, thus you must remove 
 any errors on your side. You may use the requirements.txt file to 
 install cpplint on your system.
 - docker-tests: the code is compiled within a docker, then several 
 instances are tested and the objective of the solution is checked. 
 You may create a python3 environment and install the requirements. 
 Then, you can use the script "run-actions.py" that will execute all
 the test. You may even use the option "--docker" to run the tests 
 within a docker container and emulate an environment similar to 
 github action.
 
Finally, at least the last commit must be signed in order to be merged.
To set up a signed commit, you may refer to https://docs.github.com/en/authentication/managing-commit-signature-verification
