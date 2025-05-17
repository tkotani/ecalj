import sys
from pathlib import Path
import shutil
import argparse

def main(base_dir):
    
    # POSCARALLディレクトリのパス
    poscarall_dir = Path(base_dir)        / "init0" / "POSCARALL"
    qsub_file = Path(base_dir)            / "init0" / "qsub.0"
    qsub_dependency_file = Path(base_dir) / "init0" / "qsub.dependency.0"
    
    # 必要なファイルやディレクトリが存在するか確認
    if not poscarall_dir.exists():
        print(f"Error: {poscarall_dir} does not exist.")
        sys.exit(1)
    if not qsub_file.exists():
        print(f"Error: {qsub_file} does not exist.")
        sys.exit(1)
    if not qsub_dependency_file.exists():
        print(f"Error: {qsub_dependency_file} does not exist.")
        sys.exit(1)
    
    # POSCARALL内のすべてのPOSCAR.foobarファイルを取得
    for poscar_file in poscarall_dir.glob("POSCAR.*"):
        # ファイル名から拡張子部分を取得 (foobar)
        foobar = poscar_file.name.split(".")[1]
        
        # 対応するディレクトリを作成
        target_dir = Path(base_dir) / "init" / foobar
        target_dir.mkdir(parents=True, exist_ok=True)
        
        # POSCAR.foobarファイルをコピー
        shutil.copy(poscar_file, target_dir / poscar_file.name)
        #print(f"Copied {poscar_file} to {target_dir / poscar_file.name}")
        
        # qsub.0ファイルをコピー
        shutil.copy(qsub_file, target_dir / "qsub.0")
        print(f"Copied {qsub_file} to {target_dir / 'qsub.0 ...'}")
        
        # qsub.dependency.0ファイルをコピー
        shutil.copy(qsub_dependency_file, target_dir / "qsub.dependency.0")
        #print(f"Copied {qsub_dependency_file} to {target_dir / 'qsub.dependency.0'}")

if __name__ == "__main__":
    # コマンドライン引数からベースディレクトリを取得
    parser = argparse.ArgumentParser()
    parser.add_argument("--ldir", required=True, help="LOCAL/foobar/")
    args = parser.parse_args()
    base_directory = args.ldir
    main(base_directory)
