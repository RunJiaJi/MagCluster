#!/usr/bin/env python3



def main():
    from .controller import controller
    from .args import get_magcluster_parser

    args = get_magcluster_parser()
    controller(args, subparser_name=args.subparser_name)

__all__ = [main]
#############################################
if __name__ == '__main__':
    main()