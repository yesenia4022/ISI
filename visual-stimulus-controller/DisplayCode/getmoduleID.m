function modID = getmoduleID

global GUIhandles

modID = get(GUIhandles.param.module,'value');
switch modID
    case 1
        modID = 'PG';
    case 2
        modID = 'FG';
    case 3
        modID = 'RD';
    case 4
        modID = 'FN';
    case 5
        modID = 'MP';
    case 6
        modID = 'CM';
    case 7
        modID = 'GC';
    case 8
        modID = 'SI';
    case 9
        modID = 'NM';
    case 10
        modID = 'US'; % yb: 'uncertainty stimulus'
      
end