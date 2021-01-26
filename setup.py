"""Set up files
"""

import setuptools

setuptools.setup(name="pydatadarn",
	version = "0.0.0",
	description = "Tools for analysing superdarn data",
	author = "Elliott Day",
	url = "https://github.com/ellioday/pydatadarn",
	packages = setuptools.find_packages(exclude=["test"]),
	classifiers = ["Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
	],
	python_requires=">=3.6",
	install_requires=["numpy", "scipy", "matplotlib", "spacepy", "aacgmv2"],
)
