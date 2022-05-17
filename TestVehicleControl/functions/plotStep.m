function plotStep(Map,realCarStates,k,ax)
    
    plot(ax,Map.x,Map.y,'k.-');
    hold(ax,'on');
    plotVehiclePose(ax, realCarStates(k).x, realCarStates(k).y, realCarStates(k).yaw);
    xlim(ax,[-1,5]);
    hold(ax, 'off');
    pause(0.2);

end

