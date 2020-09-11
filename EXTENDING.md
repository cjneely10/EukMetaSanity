Users with Python programming experience can easily extend **EukMetaSanity** to include additional pipelines.

Read the following explanation of how EukMS builds its workflow:

## The EukMS workflow

**EukMetaSanity** consists of four primary classes:

```
Task
TaskList
TaskManager
Data
```

`Task` provides the exact implementation needed to run a program - the specific command-line arguments needed to run a 
program.

`TaskList` is a container class that holds and calls all `Task` objects using the user-specified number of threads and 
workers.

`TaskManager` provides the order of tasks that must be done to complete a pipeline.

`Data` connects the `Task*` class functionality to the config-file with which users will interact with.

### Example:

Let's say we want to create a new pipeline that simply maps FASTQ files to a genome.