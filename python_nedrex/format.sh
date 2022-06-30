black -l 120 python_nedrex
isort --profile black python_nedrex
flake8 --max-line-length=120 python_nedrex
pylint --max-line-length=120 python_nedrex
bandit python_nedrex
mypy --strict python_nedrex
