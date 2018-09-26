lock:
	pipenv lock -r > requirements.txt
	pipenv lock -r > requirements-dev.txt
	pipenv lock --dev -r >> requirements-dev.txt


flake:
	pipenv run pyflakes .