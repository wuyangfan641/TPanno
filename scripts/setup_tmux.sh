#!/bin/bash

# 安装tmux（如果未安装）
if ! command -v tmux &> /dev/null; then
    echo "正在安装tmux..."
    if command -v apt-get &> /dev/null; then
        sudo apt-get update && sudo apt-get install -y tmux
    elif command -v yum &> /dev/null; then
        sudo yum install -y tmux
    elif command -v brew &> /dev/null; then
        brew install tmux
    else
        echo "错误：未找到支持的包管理器，请手动安装tmux。"
        exit 1
    fi
fi

# 配置别名
shell_config="$HOME/.bashrc"

# 添加别名到配置文件
cat >> "$shell_config" <<'EOF'

# Tmux 别名配置
# 1 新建会话
alias tnew="tmux new -s "
# 2 分离会话(ctrl+b d)
alias tdetach="tmux detach"
# 3 列出会话(ctrl+b s)
alias tlist="tmux ls"
# 4 接入会话
alias tattach="tmux attach -t "
# 5 杀死指定会话
alias tkill="tmux kill-session -t "
# 6 杀死全部会话
alias tkillall="tmux kill-server"
# 7 切换会话
alias tswitch="tmux switch -t "
# 8 重命名会话(ctrl+b $)
alias trename="tmux rename-session -t "
# 9 窗口上下划分窗格
alias tsplitud="tmux split-window"
# 10 窗口左右划分窗格
alias tsplitlr="tmux split-window -h"
# 11 光标到上方窗格
alias tmoveu="tmux select-pane -U"
# 12 光标到下方窗格
alias tmoved="tmux select-pane -D"
# 13 光标到左方窗格
alias tmovel="tmux select-pane -L"
# 14 光标到右方窗格
alias tmover="tmux select-pane -R"
# 15 交换窗格位置(当前窗格上移)
alias tswapu="tmux swap-pane -U"
# 16 交换窗格位置(当前窗格下移)
alias tswapd="tmux swap-pane -D"
EOF

# 提示信息
echo "Tmux配置已完成！请执行以下命令使配置生效："
echo "source $shell_config"
echo "或重新打开终端窗口。"
