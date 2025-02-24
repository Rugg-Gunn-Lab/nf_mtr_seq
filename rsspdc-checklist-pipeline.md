# Pipelines & Workflows

*[Research Software Sharing, Publication, & Distribution Checklist](https://gitlab.com/HDBI/data-management/checklists)*

Considerations for publishing a data analysis pipeline which may be used in research or as a researcher

## 📒Source control

*How can you keep track of the history of your project and collaborate on it?*

- [x] Uses [git](https://git-scm.com/) (or other source control tool such as [fossil](https://fossil-scm.org/home/doc/trunk/www/index.wiki))
	- [x]  	🥉Bronze *(Easy)*: Using version control but has a shallow project history, just placed in git for distribution
	- [ ]  	🥈Silver *(Intermediate)*: Longer project history, commit messages of mixed quality, some large messy changes
	- [ ]  	🥇Gold *(Hard)*: Silver plus - Well written commit messages, nice granular commits making discrete self-contained changes.
			Tags, releases, or branches at major project milestones, maybe some contributions from other users
	- [ ]  	🏆Platinum *(MAXIMUM OVERKILL)*: Gold plus - Some from:
			[conventional commits](https://www.conventionalcommits.org/en/v1.0.0/); Clean history with a consistent rebasing/merging strategy; Signed commits from all contributors; Contributions go through a consistent workflow like, issues, then a pull request from a branch.

<details>
<ul>
	<li>
	If the workflow tool you are using has a convenient tool for initiating a template for a workflow then you may want to get your project's git repository started using that tool.
	[Nextflow](https://www.nextflow.io/) for example has the [nf-core template](https://nf-co.re/docs/contributing/guidelines/requirements/use_the_template) which makes the creation of an nf-core style pipeline project easy.
	[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) similarly has a [standard project structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) and [template](https://github.com/snakemake-workflows/snakemake-workflow-template).
	</li>
</ul>
</details>

## ©Licencing

*On what terms can others use your code, and how can you communicate this?*

- [x] Project is suitably licensed
	- [x]  	🥉Bronze *(easy)*: There is a LICENSE file in the repository for a license which meets one of the [OSI](https://opensource.org/osd), [Debian](https://www.debian.org/intro/free), or [FSF/GNU](https://www.gnu.org/philosophy/free-sw.en.html) definitions of free/libre or open source software. Or for any contents that are not software a [Creative Commons](https://creativecommons.org/) license.
	- [ ]  	🥈Silver *(easy)*: If any prose/documentation or images is licensed differently from the code in the project this is indicated and those licences provided.
			If licences have an attribution requirement there is are easy to copy text/links for appropriate attribution.
	- [ ]  	🥇Gold *(intermediate)*: Uses [REUSE.software](https://reuse.software/) to provide license information for every file.
	- [ ]  	🏆Platinum *(intermediate)*: all previous tiers plus any images have licensing information embedded in their metadata.

<details>
<ul>
	<li>
	If you are including external code in your pipeline then you should check that their licenses are compatible and you are legally allowed to distribute your code together in this way.
	Checkout this [resource on license compatibility](https://the-turing-way.netlify.app/reproducible-research/licensing/licensing-compatibility.html).
	Generally in a pipeline your are distributing code 'alongside' other packages in a way that strong copyleft licences like the GPL intend to permit see the: [GPL FAQ](https://web.archive.org/web/20240601022910/https://www.gnu.org/licenses/gpl-faq.html#GPLInProprietarySystem)
	A pipeline that is not modifying such a library but merely using it 'as is', therefore is not considered a derivative work by these licences themselves, thus using them in this fashion is common practice, but the precise defintion of what consitutes a derived work for the purposes copyright law is generally decided on a case by case basis.
	Pipelines are sometimes run as part of the backends of web-based services so the use of an AGPL license over a GPL license may be advisable if your intent is to maximise the applicability of copy-left terms across possible use cases for your pipeline.
	</li>
	<li>
	**All software needs a license if you want to permit others to reuse it.**
	It is important to give some thought to the type of license which best suits your project, it is a choice which can have significant long term implications.
	Checkout [the turing way chapter on licensing for an introduction to the subject](https://the-turing-way.netlify.app/reproducible-research/licensing.html).
	If you have no time some pretty safe choices are: For a permissive license, the Apache 2.0.
	This would allow the re-use of your work in closed commercial code.
	For a 'copyleft' license, the GPLv3 (AGPL for server-side apps).
	This requires that anyone distributing software containing your code or derivatives of it share the source code with the people they distributed it too.
	</li>
</ul>
</details>

## 📖Documentation

*How do people know what your project is, how to use it and how to contribute?*

- [ ] Project has suitable documentation
	- [ ] 	🥉Bronze *(easy)*: Project has README file
		- [ ] README provides a description of the project structure so that the user knows which directories to find things in
	- [ ]  	🥈Silver *(intermediate)*:  If your language has a in-line documentation with which you can build documentation for your functions/objects etc. you are making full use of it.
			So that features such as the ability to produced a reference manual for all the key functional units of your code, and integrations with IDE (integrated development environment) features such as search, templates, and autocomplete can be made the most of.
	- [ ]  	🥇 Gold *(time consuming)*: Your Documentation features examples; rich, varied and working ones.
			Examples should cover both simple and more complex use cases, touching on as much of the functionality as they can.
	- [ ]  	🏆 Platinum *(intermediate)*: All examples are executed upon a full build of the documentation so that you know that they at least run.

<details>
<ul>
	<li>
		README / Manual
		<ul>
			<li>What your project is and what it does</li>
			<li>install instructions</li>
			<li>What inputs does your pipline expect and in what format?</li>
			<li>Examples of how to run the pipeline</li>
			<li>What options / flags are required / can optionally be used by your pipline?</li>
			<li>What tools does your pipeline make use of?, can it output their versions and citations</li>
			<li>What outputs does your pipline generate and in what format?</li>
			<li>
				Contribution guidance
				<ul>
					<li>for example: Issue templates, a code of conduct, process details</li>
				</ul>
			</li>
			<li>
				development environment setup
				<ul>
					<li>overview of project organisation and structure</li>
				</ul>
			</li>
		</ul>
	</li>
	<li>
	CHANGELOG it can be a good idea to include a CHANGELOG file in your project documenting things which have changed since the previous release.
	This manuscript on [The impact of package selection and versioning on single-cell RNA-seq analysis](https://doi.org/10.1101/2024.04.04.588111) provides a nice case study for why this can be useful in academic settings especially if decisions have been made to change defaults between versions.
	</li>
	<li>
		'docstrings' and similar
		<ul>
			<li>
				Many programming languages have a way of documenting your code inline which can automate the generation of some parts of the documentation.
				This often takes the form of specifically marked-up comments.
				Examples include python's [dockstrings](https://peps.python.org/pep-0257/), R's [Roxygen2](https://roxygen2.r-lib.org/), and perl's [POD](https://perldoc.perl.org/perlpod)
			</li>
		</ul>
	</li>
	<li>
		Vignettes / Examples
		<ul>
			<li>
				Examples of use of the code in the context of a real problem, beyond simple example snippets which might be included in the documentation of individual functions/objects.
				(These can also serve as a form of simple integration tests if you run them as a part of your documentation build.)
			</li>
		</ul>
	</li>
	<li>Larger projects might also include Project Documentation: Plans, Design documents and Specifications</li>
	<li>
		Process Documentation: how to proceeded with various tasks related to the project this might include: submitting issues, submitting merge requests, reporting possible vulnerabilities, testing, documentation, release, code review review
	</li>
</ul>
</details>

## 🔗Making Citable

*How should people make reference to your project and credit your work?*

- [x] Pipeline is Citable
	- [x]  	🥉Bronze *(easy)*: A [CITATION.cff](https://citation-file-format.github.io/) file exists in the code repository to provide citational metadata about your project
	- [ ]  	🥈Silver *(easy)*: bronze plus the project has persistent resolvable identifier such as a DOI, with which it can be referenced, which has been minted for the project using a tool like [zenodo](https://zenodo.org/) or [Software Heritage's Archive](https://www.softwareheritage.org/) to store an archival the project.
	- [ ]  	🥇Gold *(intermediate)*: silver plus:
		- [ ] Contributions are credited using a suitable contributor roles ontology or taxonomy (CROT) such as [CrediT](https://credit.niso.org/),[ScoRo](http://www.sparontologies.net/ontologies/scoro),[CRO](https://github.com/data2health/contributor-role-ontology), or [TaDiRAH](https://tadirah.info/).
		- [ ] All contributors are identified by their [ORCID](https://orcid.org/) or other suitable persistent identifier
	- [ ]  	🏆Platinum *(mixed)*: gold plus any two or more from:
		- [ ] (easy) All research institutions are identified by their [ROR](https://ror.org/) ID
		- [ ] (intermediate) Versioned persistent identifier with automation to update snapshots on zenodo or similar tool when a new version is created.
		- [ ] (intermediate) Annotating work cited in this work with the [Citation Typing Ontology (CiTO)](http://purl.org/spar/cito)
		- [ ] (hard) Package your software for Nix or Guix - this might not seem like it contributes to making software more citable see details below for why this is the case.

<details>

<ul>
	<li>
	Beyond merely making it possible to consistently reference a research output the higher tiers on the checklist don't merely make the output citable but help to follow better citation and bibliographic practices.
	This extends from the practical, making it easy import into reference mangers like [zotero](https://www.zotero.org/);
	To protecting against link-rot through the use of persistent digital object identifiers;
	To the use of linked / semantic data practices identify and connect contributors, the nature of their contributions, their institutional assoications and thing to which they contributed.
	Further reading on the [ethics of CROTs](https://doi.org/10.1080/08989621.2022.2161049), as well as their [evolution and adoption](https://doi.org/10.1002/leap.1496).
	</li>
	<li>
		When choosing to where to archive your code to [Zenodo](https://zenodo.org/) and/or [Software Herttage](https://www.softwareheritage.org/), they are not mutually exclusive, here is what you should know about how they differ.
		<ul>
			<li>Zenodo only stores a snapshot of your project at the time you take one and uses DOIs as persistent identifiers. </li>
			<li>Software Heritage archives your entire git history and has git commit level granularity in it's persistent identifiers [SWHIDs](https://docs.softwareheritage.org/devel/swh-model/persistent-identifiers.html) as they are partly content based i.e. derived from a hash of the contents of your repository.</li>
		</ul>
		Acording to [DataCite](https://support.datacite.org/docs/versioning) recomendations, ideally you should only take new snapshots and mint a new DOI when you create a new major version of your code, using the [semantic versioning](https://semver.org/) concepts for incrementing verions.
		They advise to only update only to update the metadata with minor and patch releases and not to mint a new DOI.
		Unfortunately the default behaviour when using the github/Zenodo integration is to mint a new DOI every time you create a new release which is typically more frequently than when you increment you major version.
		A way around this it to use a CI/CD (continious integration and deployment) job which only triggers when the major version is incremented in a release and sends a snapshot to Zenodo via their API.
		If using [gitlab2zenodo](https://gitlab.com/sbeniamine/gitlab2zenodo) achieving this would be a relatively minor modification to their suggested gitlab CI/CD workflow, one way to do this would be to store the previous version in a variable you can update update with the gitlab API and check against that.
		Snapshotting on every release could add a lot unecessarily identifiers to the DOI namespace if the project it very active and should be avoided in this case, if your project will only needs infrequent updates this option is fine.
		To maximise reproducibility it is desireable to cite software used in a given analysis with commit level version granularity, and that only versions of software that are archived be used in published research.
	</li>
	<li>
		Nix and Guix
		<ul>
			<li>
				General software repositories may not make specific provision for citation of software packages in the academic fashion.
				However some provide, what is for some use cases, a superior form of 'citation' of their own sources i.e. a complete 'software bill of materials (SBOM)'.
				This is a list of all the code used in another piece of code, its dependencies, and their dependencies recursively, along with all of their versions.
				For example [Nix](https://nixos.org/) can do this but [Guix](https://guix.gnu.org/) is perhaps the most comprehensive in its approach, though for a smaller catalog of software.
				Whilst language specific automated environment managment tools such as {[renv](https://rstudio.github.io/renv/index.html)} & [poetry](https://python-poetry.org/) provide a version of this,
				their scope does not extend to system dependencies, so whilst helpful they are incomplete.
				Guix not only provides all information necessary for a complete 'SBOM' but, it can [bootstrap](https://guix.gnu.org/en/manual/en/html_node/Bootstrapping.html) software packages in its repository from source with an extremely minimal fixed set of binaries,
				with [bitwise binary reproducibility](https://reproducible-builds.org/) for vast majority of packages, an important capability for [creating somewhat trustworthy builds](https://www.cs.cmu.edu/~rdriley/487/papers/Thompson_1984_ReflectionsonTrustingTrust.pdf).
				This creates a compute environment which is not only reproducible but 'verifiable', meaning the source of all of an environment's dependencies can in theory be scrutinised.
				It also adopts an approach to commit signing and authorisation of signers that gives it a [currently uniquely complete supply chain security architecture](https://doi.org/10.22152/programming-journal.org/2023/7/1).
				Packages or 'derivations' are 'pure functions' in the sense that only their inputs effect their outputs and they have no side-effects,
				package builds are sandboxed to prevent dependencies on any external source not explicitly provided as an input,
				and inputs are hashed to ensure that they cannot differ for the value expected when they were packaged.
				This gives these technologies an unrivaled ability to readily demonstrate the reproducibility and provenance of compute environments specified using them.
			</li>
			<li>
				In addition Guix can automatically fall back on sources archived by [software heritage](https://www.softwareheritage.org/) if the original source repository is unavailable.
			</li>
			<li>
				Example of the use of Guix in practice for bioinformatic pipelines: [PiGx](https://bioinformatics.mdc-berlin.de/pigx/) ([paper](https://doi.org/10.1093/gigascience/giy123))
			</li>
			<li>
				The primary limitation of these tools at present is not all of the software commonly used in niche scientific domains is packaged in the upstream package repositories for Nix (~100,000 packages) & Guix (~20,000 packages) despite their considerable size.
				In order to use these tools with software not yet packaged for them users must learn how to write package derivations for the tool in question and any of it's unpackaged dependencies.
			</li>
		</ul>
	</li>
</ul>
</details>

## ✅Testing

*How can you test your project so you can be confident it does what you think it does?*

- [ ] Pipeline has been tested
	- [ ]  	🥉Bronze *(easy)*: The pipeline runs and produces the expected result on a minimal example dataset
	- [ ]  	🥈Silver *(easy)*: Bronze plus:
		- [ ] A wide variety of possible input parameter combinations are tested
		- [ ] Datasets with different expected outcomes are tested, for example an output value is expected to be low in one dataset and high in another
		This provides a sense check that the pipeline is doing what you think its doing, an input perturbation produces predictable change in the output
	- [ ]  	🥇Gold *(intermediate)*: Silver Plus:
		- [ ] Tests of error handling, pipeline fails early and gracefully when given invalid inputs as parameters
		- [ ] Any Quality Control steps that you use have datasets which elicit both a pass and fail of that check
		- [ ] Test are run automatically in continuous integration tools as changes are made to the code repository
	- [ ]  	🏆Platinum *(hard)*: range of supported environments tested on different base systems with different versions in CI pipeline(s) to catch platform specific bugs

<details>
The types of testing that it may make sense to emphasise in analysis pipelines are:
Integration testing, do all the parts work together as expected? (The phrase 'expected result' is a bit nebulous, its robustness in testing hinges on how specific your expectations are.)
Some aspects might lend themselves to unit testing, but much of this may reside in the individual tools that a pipeline wraps.
Regression testing can be effective when refactoring a pipeline which is intended to produce identical results by different means.

See the software packages checklist for more details on general unit testing.

An area to focus on testing might be your quality control (QC) steps.
Checking using data that is an example of a common source data quality issues and verifying that this triggers the appropriate QC warnings, for example.
Making sure that good quality data does not trigger quality warnings.
In addition to poor quality data test: corrupt, truncated, or incorrectly formed data/configuration files.
Edge cases like many extreme values.

Test data should follow the general principle: **as small as possible as large as necessary**
Real world data can potentially be down sampled, and synthetic data generated to test edge cases.
You may be able to find a repository of pre-existing test dataset suitable for your domain nf-core provides a resource of [test datasets](https://github.com/nf-core/test-datasets), these are at least currently largely of biological data types reflecting the current composition of the Nextflow userbase.

Testing portability, does it run in different environments? nf-core nextflow pipelines for instance should be able to run with docker, apptainer (singularity), or conda based environments - do all of these work and produce the same output?

[nf-test](https://www.nf-test.com/) provides an automated testing framework for Nextflow pipelines described in [Improving the Reliability and Quality of Nextflow Pipelines with nf-test](https://doi.org/10.1101/2024.05.25.595877).
nf-test implements a number of optimisations including:
Git integration and dependeny graph analysis of pipeline component inputs/outputs to select and prioritise tests to run based on what has changed.
Convenient snapshot based regression testing.
Parallel test execution.
These features make the sometimes lengthy and laborious process of testing computationally intensive pipelines much more manageable.


[This conference paper](https://doi.org/10.1145/3322790.3330595) introduces a name for a category of software test that they call 'scientific tests' which are black-box style tests agnositic to the implementation which test the broard behaviour of a system is as predicted.
</details>

## 🤖 **Automation**

*What tasks can you automate to increase consistency and reduce manual work?*

- [ ] 
	- [ ]  	🥉Bronze *(easy)*: 
	- [ ]  	🥈Silver *(easy)*: 
	- [ ]  	🥇Gold *(intermediate)*: 
	- [ ]  	🏆Platinum *(intermediate)*: 

<details>
<ul>
	<li>
	</li>
<ul>
</details>

## 👥Peer review / Code Review

How can you get third party endorsement of and expert feedback on your project?

- [x] Pipeline has been appropriately reviewed
	- [x]  	🥉Bronze *(easy)*: Someone other than you has checked over your pipeline and given you feedback
	- [ ]  	🥈Silver *(intermediate)*: Someone other than you has successfully run your pipeline on different compute infrastructure and got the same results with test data
	- [ ]  	🥇Gold *(intermediate)*: You have published your pipeline in a pipeline repository which performs reviews of submissions such as nf-core
	  (The snakmake workflow catalog's review practices would not be adequate for this purpose)
	- [ ]  	🏆Platinum *(hard)*: Silver, gold and some from:
		- Published a peer reviewed article with a scientific review of the theoretical / statistical / mathematical underpinnings of the tool that you implemented in addition to a technical peer review of the code quality.
		  (These may well be separate reviews for example by a methods journal and a software repository reflecting their different expertise)
		- You have had and independent 'red team' attempt to find errors in your project and incorporated any relevant changes as a result
		- Your project is a part of a bug bounty program

<details>
<ul>
	<li>
		The design of the analysis, any methodological choices made and any original steps added might warrant a conventional scientific publication if for example you are making a pipeline which automates a portion of the analysis of some new datatypee
		<ul>
			<li>review the theory</li>
		</ul>
	</li>
	<li>
		Technical
		<ul>
			<li>review the implementation</li>
		</ul>
	</li>
</ul>
</details>

## 📦Distribution

How can people install or access the software emerging from your project?

- [x] Pipeline is distributed in appropriate format(s)
	- [x]  	🥉Bronze *(easy)*: Pipeline is in a software forge (such as [GitHub](https://github.com/) or [Codeberg](https://codeberg.org/)) in a standard package format so that it can be run with the pipeline manager's standard tooling.
	- [ ]  	🥈Silver *(easy)*: The software environment(s) needed by each step in your pipeline are defined so that they can be installed automatically using a suitable package / environment management tool.
	- [ ]  	🥇Gold *(intermediate)*: Pipeline is in a repository of workflows and pipelines such as [workflowhub](https://workflowhub.eu) (can be pre-release in nf-core).
	- [ ]  	🏆Platinum *(intermediate)*: Pipeline is in a curated package repository where it has undergone review and testing such as [nf-core](https://nf-co.re/).

<details>
<ul>
<li>
	Including your pipeline in a collection of pipelines increases it's visibility, can help to attract contributors, and in the case of a curated collection with good standards for how they package their pipelines provide users with confidence to they will be able to use your pipline on their compute infrastructure.
</li>
<li>
Good places to distribute workflows include:
	<ul>
	<li>
	[workflowhub](https://workflowhub.eu) most generic accepts workflows in a number different tools
	</li>
	<li>
	For Nextflow pipelines in the nf-core format [nf-core](https://nf-co.re/) (MIT license required for pipeline code)
	</li>
	<li>
	For Snakemake pipelines which conform to some relatively simple requirements [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
	</li>
	<li>
	Pipelines in the R {[targets](https://docs.ropensci.org/targets/)} tool [targetopia](https://wlandau.github.io/targetopia/)
	</li>
	</ul>
</li>
</ul>
</details>

## 💽Environment Management / Portability

How can people get specific versions of your software running on their systems?

- [x] It is possible to run your pipeline on other systems with some degree of reproducibility
	- [x]  	🥉Bronze *(easy)*: The pipeline may require some manual changes, and/or manual steps to install the requisite software environments or retrieve data, and perform some standard pre-processing of inputs like building genome references etc.
	- [ ]  	🥈Silver *(intermediate)*: External inputs can be automatically retrieved using their identifiers and pre-processing steps for these inputs are performed as apart of the pipeline.
		At least one method for describing the required compute environment(s) is supplied e.g. conda environments, or container build files so that others can be specified if a user has compute infrastructure which does not support the method(s) supplied.
	- [ ]  	🥇Gold *(hard)*: Complete pipeline can be executed with single command with all data and dependencies fetched automatically, barring the need for any system specific resource constraint configuration.
	- [ ]  	🏆Platinum *(practically impossible in practice)*: The entire dependency tree of your pipeline including the pipeline manager itself, the OS it's running on and the firmware of the hardware it's running on can be bootstraped from source and produce bitwise identical binaries, also including any pre-processing of data inputs to the pipeline such as base calling for sequencing data.

<details>
<ul>
<li>
In the context of a pipeline each independent step should ideally be performed in it's own environment , perhaps defined within a container, with only the tools necessary to perform that step of the analysis.
Many pipeline management tools support specifying per-task compute environments using tools such as conda, and container technologies such as docker and singularity/apptainer.
	<ul>
		<li>
			<p>
			Container images are a convenient format in which to distribute software along with its dependencies and to isolate this environment from other software which can help avoid any conflicts in dependencies.
			Container images like virtual machine images can be quite large and thus cumbersome to distribute, they are also something of a black box once built.
			Containers provide many of the advantages of virtual machines (VM) but generally with less performance penalty.
			However whilst they are narrowly reproducible they are not readily interrogated and checked unless you provide the build instructions which generated the image, for example a Dockerfile.
			Unfortunately the process of building container images is itself not necessarily reproducible.
			Thus when specifying container builds it is best practice to specify exact package versions in your build so that the image builds are reproducible.
			This can be challenging as many popular operating system package managers lack the tooling to do this easily.
			Pinning your container build to a snapshot of package repositories taken at a given date that will be available archivally is one way to address this.
			This paper provides some [Recommendations for the packaging and containerizing of bioinformatics software](https:/doi.org/10.12688/f1000research.15140.2)
			In the case of bioinformatics pipelines it is often easiest to specify your environment with [conda](https://conda.org/) and then build containers and/or VMs which install that same conda environment on linux base image such as debian.
			This has the advantage that the conda environment can be used independent of any images built with it reducing the maintenance burden for supporting multiple approaches to distributing the compute environment.
			Which package/environment management tool has well packaged versions of all the relevant software may be specific to your discipline.
			</p>
			<p>
			If you do need to inspect the contents of a container images a number of tools developed by [Anchore](https://github.com/anchore) in particular [Syft](https://github.com/anchore/syft) can be helpful in producing an account of the software installed in the image.
			It is a best practice to keep the software that you install in an image to the minimum necessary for the function that you need the image to perform, however, determining what this minimal set is can be non-trivial.
			</p>
		</li>
		<li>
			<p>
			Functional package managers such as [Nix](https://nixos.org/) and [Guix](https://guix.gnu.org/) have a 'best practices by design' approach to packaging software.
			They do not suffer from the issue of it being difficult to determine what is and is not a required dependency as this work is done up-front when the software is packaged.
			They usually require that dependencies be completely specified and packages be build in a sandboxed environment which only has access to the explicitly specified dependencies.
			This provides much stronger guarantees of the ability to specify and build reproducible environments.
			It is also possible to build container and VM images specified with these tools, and a container specified with them could be a drop in replacement for one specified with conda and Docker for example.
			Unfortunately these tools have yet to see wide adoption in the scientific / research computing communities and thus many packages used by these communities are not packaged for these tools, (despite nixpkgs being the largest extant software package respoitory with >100,000 packages), hindering their broarder adoption.
			Use and awareness is growing and there are some excellent [case studies](https://bioinformatics.mdc-berlin.de/pigx/).
			Nix is also cross platform working natively on MacOS, on windows via the windows subsystem for linux and even on android.
			</p>
			<p>
			It is worth being aware of these tools and considering packaging any software that you produce for them as they are gaining popularity and address many of the shortcommings and limitations of current package and environment mangement solutions.
			</p>
		</li>
	</ul>
</li>
	<li>
	Are the environments for each step of your pipeline well described using an environment management tool such as [Conda](https://conda.org/), or [Spack](https://spack.readthedocs.io/en/latest/) and/or supplied as [OCI](https://opencontainers.org/) containers, runnable with tools such [Docker](https://www.docker.com/), [podman](https://podman.io/), [lxc](https://linuxcontainers.org/lxc/introduction/), [Singularity/Apptainer](https://apptainer.org/), or others?
	</li>
	<li>
	Many popular pipeline management tools integrate with environment management and container runtimes to facilitate portability of reproducible compute environments. see:
	</li>
	<ul>
		<li>
		[Snakemake - integrated Package management](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management)
		</li>
		<li>
		[Nextflow - containers](https://www.nextflow.io/docs/latest/container.html)
		</li>
	</ul>
	<li>
	The Pipelines in Genomics [PiGx](http://bioinformatics.mdc-berlin.de/pigx/) ([paper](https://doi.org/10.1093/gigascience/giy123)) collection represents a gold standard in reproducible computational environments for genomics pipelines.
	It used the [Guix](https://guix.gnu.org/) functional package manager to attain >97% bitwise reproducibility for dependencies across the pipelines in their collection.
	Containers built with Nix or Guix can be used in pipeline managers, as PiGx does with snakemake.
	</li>
</ul>
</details>

## 🌱 Energy Efficiency

How can you and your users minimise wasted energy?

- [ ] Consideration has been given to the energy efficiency of the pipeline
	- [ ]  	🥉Bronze *(intermediate)*:
		- [ ] Adopt a process to review outputs and keep only those necessary
	- [ ]  	🥈Silver *(intermediate)*: Bronze plus:
		- [ ] Energy consumption and emissions Estimation/Measurement
	- [ ]  	🥇Gold (hard) silver plus, some from:
		- [ ] Appropriately balances the possible efficiency gains against the cost of testing, development, and benchmarking
		- [ ] Makes use of efficient algorithms for its task
		- [ ] Makes use of tools implemented in an efficient language for the task
	- [ ]  	🏆Platinum (hard) gold plus, some from:
		- [ ] Makes use of hardware acceleration if appropriate
		- [ ] Can make use of carbon intensity aware run scheduling
		- [ ] Binaries are compiled with flags optimal for performance on the architecture/hardware on which the computation is performed if a job is to be repeated enough times to make sense given the energy expended on compilation.

<details>

Everyone likes fast and efficient code, but especially if your code is going to be re-used by a lots of people in a computationally demanding application it can consume a lot of energy.
This translates to carbon emissions, water use and opportunity costs for whatever else could have been done with that energy and compute time.

Consider what can you do to make your code a little more efficient:

<ul>
	<li>
	Measurement - make use of tools to estimate the energy &/or carbon costs of your analyses and report them to the end user.
	If you are using a [Nextflow](https://www.nextflow.io/) pipeline then there is a very easy way to estimate the energy utilisation and carbon footprint of your pipelines with [`nf-co2footprint`](https://nextflow-io.github.io/nf-co2footprint/).
	This is a generalisation of a tool originally developed to estimate emissions for jobs submitted to HPC clusters running the [SLURM](https://slurm.schedmd.com/overview.html) scheduler: [GA4HPC](https://www.green-algorithms.org/GA4HPC/), see: [Green Algorithms: Quantifying the Carbon Footprint of Computation](https://doi.org/10.1002/advs.202100707).
	The people in the [Green Algorithms](https://www.green-algorithms.org/) community of practice have some useful advice and resources for anyone interested in this subject.
	</li>
	<li>
	Good documentation and good error handing/messages can reduce the number of times people make mistakes using your code that means they re-run or partially re-run it fewer times before they figure out how to use it correctly.
	</li>
	<li>
	Don't generate unnecessary outputs that will sit on people's drives unused, clean results of intermediate steps.
	If you're making a pipeline produces a lot of intermediate files and outputs consider which of these are needed or good defaults, which could be optional and which could be discarded by default.
	Defaults are king and people will mostly keep whatever your tool outputs often essentially indefinitely so you can reduce the energy expended on unnecessary storage by keeping your outputs lean.
	You might have varying degrees of verbosity of output with a more verbose mode for debugging but defaulting to just the essentials.
	</li>
	<li>
	For pipelines in particular caching results and avoiding needing to re-compute things if possible - make best use of these features in pipeline managers for example by having small granular tasks to minimise repeated work on run failure.
	</li>
	<li>
	Choice of libraries and frameworks - some libraries may be more efficient that others or be a wrapper around an efficient implementation in another language, or be able to make use of offload to hardware accelerators.
	</li>
	<li>
	Offload to hardware accelerators where available, vector matrix and array arithmetic can often benefit from very substantial speed-ups on hardware specialised for these types of calculations, or even binaries compiled with the right instruction set extensions enabled to take full advantage of hardware acceleration features on many CPUs.
	Doing this directly can be quite challenging but using libraries capable of managing this offload for you can make it more approachable. (This can potentially introduce interesting reproducible computation challenges due to things like differences in handling of floating point arithmetic between hardware/firmware implementations.)
	</li>
	<li>
	Benchmarking & Profiling to locate and improve inefficient code.
	Don't optimise prematurely - it is often surprising which pieces of your code turnout to be slow, measure it first and check where to focus your attention.
	This can go hand in hand with having done robust testing as a good test suite means that you can confidently re-factor an inefficient piece of code without fear of introducing errors.
	Robust testing, especially of large computationally intensive pipelines can be energy intensive so it is important to make efficient test suites which still provide robust coverage.
	See the testing section for more information on optimisations in the [nf-test](https://www.nf-test.com/) framework which help to make testing Nextflow pipelines more efficient.
	</li>
	<li>
	Language Choice - some languages are (at least on average) more efficient than others, weigh this as a factor when selecting a language.
	In the context of a pipeline manger the efficiency of the language it itself is written in is rarely particularly relevant, what matters more are features which increase the efficiency of overall execution such as good management of caching and degrees of parallelism.
	However the efficiency of the individual components of the pipeline may be of disproportionate impact as they may be run many time by many people so selecting the most efficient modules is of higher impact.
	</li>
	<li>
	Consider timing and location of computationally intensive runs, you could automate running at optimal times and in optimal places to make use of surplus renewable energy. In the UK you could query the [national grid carbon intensity API](https://www.carbonintensity.org.uk/) to pick opportune times and/or places.
	</li>
</ul>
</details>


## ⚖ Governance, Conduct, & Continuity

How can you be excellent to each other, make good decisions well, and continue to do so?

- [ ] The project has a suitable governance model
	- [ ]  	🥉Bronze *(easy)*: The governance model is clearly communicated
	- [ ]  	🥈Silver *(easy)*: Bronze Plus - Project has continuity planning in place (2 or more from)
		- [ ] Source archived and/or mirrored to other platforms
		- [ ] Public archives of key project governance documentation and plans for continuity of operations in the events such as the loss of key project infrastructure
		- [ ] Plans of action in the event project admin(s) are no longer available
	- [ ]  	🥇Gold *(intermediate)*: Project has a governance model appropriate to its scale
		- [ ]
		- [ ] Project has clear and transparent processes
	- [ ]  	🏆Platinum *(Hard)*: Project has a track record of good governance and policy
		- [ ]
		- [ ] Learned from any mistakes and implemented policy changes as a result

<details>
<ul>
	<li>
	If you are the Benevolent Dictator For Life (BDFL) of your project and the Code of Conduct (CoC) is "Don't be a Dick" that's fine, for many individual hobby projects this a functional reality.
	Becoming a BDFL tends to be the default unless you take steps to avoid it and cultivate community governance as your project begins to grow - failing to do this and being stuck in charge can become quite the burden in sucessful projects.
	Be warning of adopting policies that you lack resources, time, interest, skill, or inclination to be an active enforcer, mediator and moderator of community norms and disputes, It is helpful to be clear about what you can and cannot commit to doing.
	Only by communicating this might you be able to find community members to help you with setting and enforcing these norms, if or when your community attains a scale where this becomes relevant - community management is its own skill set.
	If you can't moderate them avoid creating and/or continuing ungoverned community spaces that can become a liability for you and your project's reputation.
	Just as there are off-the-shelf licenses there are off-the-shelf codes of conduct, the [Contributor Covenant](https://www.contributor-covenant.org/) is perhaps the best known and most widely used, though may need some customisation to your needs.
	Adopting such a CoC gives you some guidance to follow if there is bad behaviour in your project's community and communicates that you as the project leadership take the responsibility of creating a respectful environment for collaboration seriously.
	It can also signal that your project is a place where everyone is welcome but expected to treat one another with respect, and that failing to do so will result in penalties potentially including exclusion from the community.
	[The Turing Way](https://doi.org/10.5281/zenodo.3233853) provides quite a nice example of a [CoC developed specifically for their project](https://web.archive.org/web/20240412122958/https://book.the-turing-way.org/)
	You will need to provide contact information for the person(s) responsible for the enforcement of the CoC in the appropriate place and be able to follow up in the event it is used.
	git forges often recognise files with the name `CODE_OF_CONDUCT.md` in the root of project and provide a link to them on project home pages, so this is a good place to document such policies.
	If you are the BDFL of a small project then interpretation and enforcement of such a CoC tends to fall solely on you - game out some courses of action for what you'd do if faced with some common moderation challenges.
		<ul>
			<li>
				Once a project attracts a larger community there is greater scope for disputes and therefore for the need for dispute resolution mechanisms.
				Free/Libre and Open Source Software development and maintenance can be thought of as a commons so I would refer you to the [work of Elinor Ostrom](https://archive.org/details/ElinorOstromGoverningTheCommons) on how commons have been successfully (or unsuccessfully) governed when thinking about what processes to adopt for your project.
				More recently [Nathan Schneider's Governable Spaces: Democratic Design for Online Life](https://doi.org/10.1525/luminos.181) tackles some of these issues as applied to online spaces.
			</li>
			<li>
				This is summarised in the [8 Principles for Managing a Commons](https://web.archive.org/web/20240411163813/https://onthecommons.org/magazine/elinor-ostroms-8-principles-managing-commmons/index.html)
				<ol>
					<li>Define clear group boundaries.</li>
					<li>Match rules governing use of common goods to local needs and conditions.</li>
					<li>Ensure that those affected by the rules can participate in modifying the rules.</li>
					<li>Make sure the rule-making rights of community members are respected by outside authorities.</li>
					<li>Develop a system, carried out by community members, for monitoring members’ behaviour.</li>
					<li>Use graduated sanctions for rule violators.</li>
					<li>Provide accessible, low-cost means for dispute resolution.</li>
					<li>Build responsibility for governing the common resource in nested tiers from the lowest level up to the entire interconnected system.</li>
				</ol>
			</li>
			<li>
				An informal [do-ocracy](https://web.archive.org/web/20230201164442/https://communityrule.info/modules/do-ocracy/) in the fiefdom of BDFL is often the default state of projects that have not given much conscious thought to how they want to be governed and are thus often subject to many of the same common failure modes of this model.
				How are decisions made in your project? Do you need the mechanisms of governance used by community and civil society organisations?
				By-laws, a committee and/or working groups, general meetings, votes, minutes? A version of these may be necessary to avoid [The Tyranny of Structurelessness](https://web.archive.org/web/20240406094310/https://www.jofreeman.com/joreen/tyranny.htm)
				How can you map these onto your development infrastructure and make the decisions of your governing bodies enactable and enforceable?
			</li>
		</ul>
	</li>
	<li>
	Continuity planning: What happens to your project if something happens to you?
	The code will likely live on due the distributed nature of git but what about the issue tracker, the website etc.
	Who else has the highest level of privilege on your project or a mechanism to attain it?
	The principle of least privilege dictates that you keep the number of people with this level of access to a minimum but you may then create a single point of failure.
	Password managers like [bitwarden](https://bitwarden.com/) have a feature where designated people can be given access to your vault if they request it and you do not deny it within a certain time-frame.
	This could provide a lower level admin with a mechanism to escalate their privileges if you are unable to do this for them.
	However, this delay might be an issue for continuity of operations if administrator action is needed within the waiting period.
	Game it out, have a plan, write it down, let people know you have a plan.
	</li>
	<li>
	Does your project take donations?
	Does it have a trademark?
	Does it need a legal entity to hold these?
	Who is on the paperwork and who has signing authority?
	Who keeps track of expenditures?
	Tools & Organisations like [OpenCollective](https://opencollective.com/) can help with some of these issues.
	</li>
	<li>
	If your project has potential cybersecurity implications what procedures do you have in place for people to disclose vulnerabilities in the project so that they can be patched before they are made public.
	What systems do you have in place to disclose a vulnerability once it has been patched and ensure that users know that they need to update.
	</li>
	<li>Whole project data longevity - what plans do you have in place to backup and archive materials pertaining to your project that are not under source control?</li>
	<li>
	User support
		<ul>
		<li>What support can users expect, or not expect?</li>
		<li>Where can they ask for it?</li>
		<li>Is there somewhere where users can provide support to other members of the user community, such as a forum?</li>
		<li>Can they pay for more support?</li>
		<ul>
	</li>
</ul>

</details>

 <p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://gitlab.com/HDBI/data-management/checklists">Research Software Sharing, Publication, & Distribution Checklists</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://orcid.org/0000-0002-2574-9611">Richard J. Acton</a> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>

