import networkx as nx
from collections import namedtuple
from typing import List, Type, Tuple, Dict
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.tasks.dependencies import dependencies

Node = namedtuple("Node", ("name", "scope", "dependency_input"))


class DependencyGraphGenerationError(BaseException):
    """ Exception class if unable to form directed acyclic (dependency) graph from
    provided pipeline configuration

    """
    pass


class DependencyGraph:
    """ Class takes list of tasks and creates dependency graph of data
    Topological sort outputs order in which tasks can be completed

    All tasks are assumed to require (at least) a root task to complete

    The root task is simply populating input files to run in a pipeline

    """

    def __init__(self, tasks: List[TaskList]):
        self.idx: Dict[str, TaskList] = {task.name: task for task in tasks}
        self.idx.update(dependencies)
        self.graph = nx.DiGraph()
        self._build_dependency_graph(tasks)
        if not nx.is_directed_acyclic_graph(self.graph):
            raise DependencyGraphGenerationError(
                "There was an error generating dependency graph information for this pipeline.\n"
                "If you are a user seeing this, this is likely a fatal error - contact the authors of the pipeline\n"
                "If you are a developer, double-check dependency information in your pipeline"
            )

    def _build_dependency_graph(self, tasks: List[TaskList]):
        root = Node(name="root", scope="", dependency_input=("root", "fna"))
        self.graph.add_node(root)
        # Add dependencies stored in TaskList object's .requires member
        for task_list in tasks:
            task_node = Node(name=task_list.name, scope="", dependency_input=("root", "fna"))
            self.graph.add_edge(root, task_node)
            for requirement in task_list.requires:
                self.graph.add_edge(Node(name=requirement, scope="", dependency_input=("root", "fna")), task_node)
            for dependency in task_list.depends:
                self._add_requirements_within_dependencies(self.graph, Node(name=dependency.name, scope=task_list.name,
                                                                            dependency_input=dependency.input),
                                                           task_node, "requires", "")

    def _add_requirements_within_dependencies(self, graph: nx.DiGraph, node: Node, task_node: Node, attr: str,
                                              scope: str):
        for requirement in getattr(self.idx[node.name], attr):
            if attr == "depends":
                new_node = Node(name=requirement.name, scope=scope, dependency_input=node.dependency_input)
            else:
                new_node = Node(name=requirement, scope=scope, dependency_input=node.dependency_input)
            graph.add_edge(new_node, task_node)
            self._add_requirements_within_dependencies(graph, new_node, task_node, attr, scope)

    def _get_dependencies_at_level(self, task: TaskList) -> List[Tuple[Type[TaskList], str, Tuple[str, str]]]:
        graph = nx.DiGraph()
        task_node = Node(name=task.name, scope="", dependency_input=("root", "fna"))
        graph.add_node(task_node)
        for dependency in task.depends:
            graph.add_edge(Node(name=dependency.name, scope=task.name, dependency_input=dependency.input), task_node)
            self._add_requirements_within_dependencies(graph, Node(name=dependency.name, scope=task.name,
                                                                   dependency_input=dependency.input), task_node,
                                                       "depends", task.name)
        sorted_nodes = list(nx.topological_sort(graph))
        sorted_nodes.remove(task_node)
        return [(self.idx[node.name], str(node.scope), node.dependency_input) for node in sorted_nodes]

    @property
    def sorted_tasks(self) -> List[Tuple[Type[TaskList], str, Tuple[str, str]]]:
        """ Run topological sort of all tasks in pipeline and output in order that allows
        for completion of dependencies in required order

        :return: List of TaskList child classes to run
        """
        sorted_nodes: List[Node] = list(nx.topological_sort(self.graph))
        sorted_nodes.remove(Node(name="root", scope="", dependency_input=("root", "fna")))
        final_sort: List[Tuple[Type[TaskList], str, Tuple[str, str]]] = []
        for node in sorted_nodes:
            final_sort += self._get_dependencies_at_level(self.idx[node.name])
            final_sort.append((self.idx[node.name], node.scope, node.dependency_input))
        return final_sort
