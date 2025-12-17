#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
SCRIPT DE BUILD PARA EXECUTÁVEL CAR-T REORDER v2.0
================================================================================

Este script cria um executável standalone do Pipeline CAR-T usando PyInstaller.

REQUISITOS:
    pip install pyinstaller

USO:
    python build_executavel.py

O executável será gerado na pasta 'dist/'
================================================================================
"""

import subprocess
import sys
import os
import shutil


def verificar_pyinstaller():
    """Verifica se o PyInstaller está instalado."""
    try:
        import PyInstaller
        print(f"PyInstaller encontrado: versão {PyInstaller.__version__}")
        return True
    except ImportError:
        print("PyInstaller não encontrado!")
        print("\nInstale com: pip install pyinstaller")
        return False


def limpar_builds_anteriores():
    """Remove arquivos de builds anteriores."""
    diretorios = ['build', 'dist', '__pycache__']
    arquivos = ['reordenacao-cart-v2-gui.spec']

    for d in diretorios:
        if os.path.exists(d):
            print(f"Removendo diretório: {d}")
            shutil.rmtree(d)

    for f in arquivos:
        if os.path.exists(f):
            print(f"Removendo arquivo: {f}")
            os.remove(f)


def criar_executavel():
    """Cria o executável usando PyInstaller."""
    print("\n" + "="*60)
    print("CRIANDO EXECUTÁVEL CAR-T REORDER v2.0")
    print("="*60)

    # Nome do arquivo principal
    arquivo_principal = "reordenacao-cart-v2-gui.py"

    if not os.path.exists(arquivo_principal):
        print(f"ERRO: Arquivo '{arquivo_principal}' não encontrado!")
        return False

    # Opções do PyInstaller
    opcoes = [
        'pyinstaller',
        '--onefile',              # Criar um único arquivo executável
        '--windowed',             # Não mostrar console (Windows)
        '--name=CAR-T_Reorder',   # Nome do executável
        '--clean',                # Limpar cache antes de compilar
        arquivo_principal
    ]

    # Adicionar ícone se existir
    if os.path.exists('icone.ico'):
        opcoes.insert(-1, '--icon=icone.ico')

    print("\nExecutando PyInstaller...")
    print(f"Comando: {' '.join(opcoes)}")
    print("\n" + "-"*60)

    try:
        resultado = subprocess.run(opcoes, check=True)
        print("-"*60)
        print("\nExecutável criado com sucesso!")

        # Verificar onde o executável foi criado
        if os.name == 'nt':  # Windows
            exe_path = os.path.join('dist', 'CAR-T_Reorder.exe')
        else:  # Linux/Mac
            exe_path = os.path.join('dist', 'CAR-T_Reorder')

        if os.path.exists(exe_path):
            tamanho = os.path.getsize(exe_path) / (1024 * 1024)  # MB
            print(f"\nLocalização: {os.path.abspath(exe_path)}")
            print(f"Tamanho: {tamanho:.1f} MB")
        else:
            print(f"\nExecutável esperado em: {exe_path}")

        return True

    except subprocess.CalledProcessError as e:
        print(f"\nERRO ao criar executável: {e}")
        return False
    except FileNotFoundError:
        print("\nERRO: PyInstaller não encontrado no PATH!")
        print("Tente: pip install pyinstaller")
        return False


def main():
    """Função principal."""
    print("="*60)
    print("BUILD SCRIPT - CAR-T REORDER v2.0")
    print("="*60)

    # Mudar para o diretório do script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"\nDiretório de trabalho: {script_dir}")

    # Verificar PyInstaller
    if not verificar_pyinstaller():
        print("\n" + "="*60)
        print("INSTALAÇÃO DO PYINSTALLER:")
        print("="*60)
        print("\n1. Execute: pip install pyinstaller")
        print("2. Execute novamente este script")
        print("\nOu execute diretamente:")
        print("   pyinstaller --onefile --windowed --name=CAR-T_Reorder reordenacao-cart-v2-gui.py")
        return 1

    # Perguntar se deseja limpar builds anteriores
    print("\n" + "-"*60)
    print("Deseja limpar builds anteriores? (s/n): ", end="")

    try:
        resposta = input().strip().lower()
        if resposta in ['s', 'sim', 'y', 'yes']:
            limpar_builds_anteriores()
    except:
        # Em caso de execução não-interativa, limpar automaticamente
        limpar_builds_anteriores()

    # Criar executável
    if criar_executavel():
        print("\n" + "="*60)
        print("BUILD CONCLUÍDO COM SUCESSO!")
        print("="*60)
        print("\nO executável está em: dist/CAR-T_Reorder")
        print("\nVocê pode copiar este arquivo para qualquer computador")
        print("e executá-lo sem precisar do Python instalado.")
        return 0
    else:
        print("\n" + "="*60)
        print("BUILD FALHOU!")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
