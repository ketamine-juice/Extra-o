# Tutorial Rápido de Git (para Windows)

[Instalar o Git](https://github.com/git-for-windows/git/releases/download/v2.44.0.windows.1/Git-2.44.0-64-bit.exe)

[Instalar o GitHub CLI - para clonar e fazer pull requests diretamente da conta](https://github.com/cli/cli/releases/download/v2.45.0/gh_2.45.0_windows_amd64.msi)

*Todos os passos assumem estar a correr o Terminal na pasta onde é pretendido que a ação corra*
*Exemplo: Navegar até à pasta onde é pretendido clonar, right-click e "Abrir Terminal Aqui"*

## Clonar um repositório

> Correr o comando `git clone https://github.com/url-do/repositorio.git` (e.g. `git clone https://github-com/ketamine-juice/Extra-o`

## Fazer pull e sync

> `git fetch`
> `git pull`

## Criar um branch novo

Criar um branch é importante para evitar commits concorrentes e quebrar ficheiros quando estamos a trabalhar em paralelo.
Evita, por exemplo, Jupyter notebooks quebrados. 

> Correr o comando `git checkout -b nome-do-branch`  

## Fazer um commit

> `git commit -m "Mensagem de commit"` 

## Fazer push das mudanças

Sem branch:
> `git push` 

Com branch:
> `git push -u origin nome-do-branch`
 
## Fazer um pull request

Quando estamos a trabalhar num branch à parte é necessário fazer um pull request para que as alterações integrem o ramo principal do repositório. Isto é feito diretamente no Github, navegando até ao repositório > pull requests > new pull request.
