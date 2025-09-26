## Workflows directory
This directory provides convinient workflows to process multiple of samples and append individual results into overall workflows. 

The `default` directory provides workflows that generate final report with the header. These are recommended workflows designed for users with admin Galaxy rights

The `usegalaxy` directory contains workflows that are designed to run on [`usegalaxy.eu`](https://usegalaxy.eu/) server respecting the context of the server. Currently reports are **without** headers

### Categories
* `Contig read mode` workflows designed for both forward and reverse read input per sample
* `Forward read mode` workflows designed for forward only reads
* `Reverse read mode` workflows designed for reverse only reads

### Installation
Workflows should be imported into individual Galaxy instances and accounts via a file import. If any of the dependencies are missing, you will be notified. Admin of the server could easily install required dependencies depending on the workflow. 

### Dependencies
The above workflows depend on one of the following dependencies to collate individual reports into a single master one. Check each workflow to see which one applies.

* `concatenate_multiple_datasets`
	* used by the `default` workflows
	* description: Concatenate multiple datasets tail-to-head, including collection datasets. 
	* owner: artbio
	* version: 1.4.1
	* revision: 7 (2019-07-09)
	* [https://toolshed.g2.bx.psu.edu/view/artbio/concatenate_multiple_datasets](https://toolshed.g2.bx.psu.edu/view/artbio/concatenate_multiple_datasets)	

* `concat_text_files`
	* description: Concatenates any text based files you select from your history
	* owner: mandorodriguez
	* version: 1.0.0
	* revision: 0 (2015-09-28) 	
	* [https://toolshed.g2.bx.psu.edu/view/mandorodriguez/concat_text_files](https://toolshed.g2.bx.psu.edu/view/mandorodriguez/concat_text_files)

* `text_processing`
	*  used by the Galaxy Europe workflows
	*  description: Concatenate datasets
	*  owner: bgruening
	*  version: 0.1.0 (2019-04-03)
	*  revision: 13
	*  [https://toolshed.g2.bx.psu.edu/repository/view_tool_metadata?repository_id=2593fd36ae8011aa&changeset_revision=0a8c6b61f0f4&tool_id=tp_cat&render_repository_actions_for=tool_shed](https://toolshed.g2.bx.psu.edu/repository/view_tool_metadata?repository_id=2593fd36ae8011aa&changeset_revision=0a8c6b61f0f4&tool_id=tp_cat&render_repository_actions_for=tool_shed)

