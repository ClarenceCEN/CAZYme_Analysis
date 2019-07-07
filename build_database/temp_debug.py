import re

def find_value_recursive(cazy_dict,cazyme):
    name_pattern = re.compile(r'(GH|GT|CBM|PL|AA|CE)(\d+)(_\d+)?')
    try:
        cazy_cat = re.search(name_pattern, cazyme).group(1)
        cazy_num = re.search(name_pattern, cazyme).group(2)
    except:
        print('Find %s!' % cazyme)
        print("Interesting family! What is that?")

    global flag
    global levels
    flag = 0
    levels = cazyme
    for key,value in cazy_dict.items():
        if isinstance(value,dict):
            global n
            n += 1
            find_value_recursive(value, cazyme)
            if flag == 1:
                n -= 1
                levels =  'L'+str(n)+'_'+key+';'+levels
                return levels
        elif cazyme == value:
            n += 1
            flag = 1
            levels =  'L'+str(n)+'_'+key+';'+'L'+str(n+1)+'_'+levels
            return





def main():

    CAZYme_family = {'CAZyme':{'Enzyme_Classes':{'GlycosylTransferases':'GT',
              'GlycosideHydrolases':'GH',
              'CarbohydrateEsterases':'CE',
              'PolysaccharideLyases':'PL',
              'AuxiliaryActivities':'AA'},
              'Associated_Modules':{'Carbohydrate-BindingModules':'CBM'}}}

    global n
    n = 0

    a  = find_value_recursive(CAZYme_family,'GH3_13')
    print(a)

if __name__  == '__main__':
    main()