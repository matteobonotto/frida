


test:
	poetry run pytest -m 'not dev' -vs

style:
	poetry run black src/frida --config pyproject.toml

type:
	poetry run mypy src/frida --config pyproject.toml