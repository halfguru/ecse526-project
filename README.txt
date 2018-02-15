- All the source code is located in the "code" folder.
- You can run the algorithm from design_controller.m.
- There is extensive commenting across all the source code.
- To change the evolution parameters, refer to design_controller.m.
- To change or add trajectories, refer to trajectory_definition.m.

Note: to prevent accidentaly closing the figures while NEAT is running, you must use the following code snippet in order to force the figures to close if you want to quit MATLAB:
figAllHandle = findall(0,'Type','figure');
if ~isempty(figAllHandle)
    for iFigHandle = 1:numel(figAllHandle)
        set(0, 'currentfigure', figAllHandle(iFigHandle));
        close('force');
        delete(figAllHandle(iFigHandle));
    end
    clear figAllHandle;
end