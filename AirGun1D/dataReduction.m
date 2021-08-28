% Reduce futures object to wall states, shuttle and bubble
% Use this function if machine has insufficient memory
function structReduced = dataReduction(solnStruct)
    structReduced.qWall = solnStruct.q(1:3,:);
    structReduced.qOut = solnStruct.q(end-2:end,:);
    structReduced.shuttle = solnStruct.shuttle;
    structReduced.bubbleContinuationTime = ...
        solnStruct.bubbleContinuationTime;
    structReduced.t = solnStruct.soln.x;
    structReduced.bubbleContinuationState = ...
        solnStruct.bubbleContinuationState;
end