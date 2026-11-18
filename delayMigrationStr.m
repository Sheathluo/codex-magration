function str = delayMigrationStr(x)
    switch x
        case 1
            str = '符号级徙动';
        case 2
            str = '采样点级徙动';
        case 3
            str = '不徙动';
        otherwise
            str = '未知模式';
    end
end