PIP=pip3

.PHONY: docs

lock:
	pipenv lock
	pipenv lock -r > requirements.txt
	pipenv lock -r > requirements-dev.txt
	pipenv lock --dev -r >> requirements-dev.txt

pylint:
	pipenv run pylint -E jdna

flake:
	pipenv run pyflakes .

docs:
	@echo "Updating docs"

	# copy README.md to README.rst format for Sphinx documentation
	# we can comment this out if we do not want to include the README.md in the sphinx documentation
	#pipenv run pandoc --from=markdown --to=rst --output=docsrc/README.rst README.md

	pandoc --from=markdown --to=rst --output=README.rst README.md
	rm -rf docs
	cd docsrc && pipenv run make html
	find docs -type f -exec chmod 444 {} \;
	@echo "\033[95m\n\nBuild successful! View the docs homepage at docs/html/index.html.\n\033[0m"

	touch docs/.nojekyll

klocs:
	find . -name '*.py' | xargs wc -l

benchmark:
	rm -rf .benchmarks/images/*svg
	pipenv run python -m pytest --benchmark-autosave --benchmark-max-time=0.1 --benchmark-group-by=func --benchmark-histogram=.benchmarks/images/histogram

