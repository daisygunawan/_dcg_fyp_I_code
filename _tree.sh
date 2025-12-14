tree -h --du -L 4 \
  -I "_devfiles|__pycache__|*.pyc|venv_proto|.git|.DS_Store" \
  --dirsfirst . > _tree_prototype.txt