## ECSE 526 final project

Control of a Quadrotor Using Generated ANN With NEAT Algorithm implement using Matlab and Simulink.

The aim of this research is control the nonlinear system  (UAV drone quadcopter commercial model) with a genetic algorithm implementation called NEAT (NeuroEvolution of Augmenting Topologies) by Stanley and Miikkulainen. The reasoning is that in general, due to their random nature, evolutionary algorithms can’t find an optimal solution but will often find a good or sufficient solution. This project seeks to observe the relevance of outputs obtained from neural network controllers. It’s worth mentioning that it is a highly experimental approach.

## Notes

* All the source code is located in the "code" folder.
* You can run the algorithm from design_controller.m.
* There is extensive commenting across all the source code.
* To change the evolution parameters, refer to design_controller.m.
* To change or add trajectories, refer to trajectory_definition.m.

To prevent accidentaly closing the figures while NEAT is running, you must use the following code snippet in order to force the figures to close if you want to quit MATLAB:
figAllHandle = findall(0,'Type','figure');
if ~isempty(figAllHandle)
    for iFigHandle = 1:numel(figAllHandle)
        set(0, 'currentfigure', figAllHandle(iFigHandle));
        close('force');
        delete(figAllHandle(iFigHandle));
    end
    clear figAllHandle;
end

![Matlab diagram](https://i.imgur.com/jf99Jat.png)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgements
* Maxence Boutet, teammate

